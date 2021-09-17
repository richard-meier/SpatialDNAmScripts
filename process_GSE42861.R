
options(stringsAsFactors=FALSE)
options(timeout = 300)
library(Biobase)
library(GEOquery)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(quadprog)
library(R2OpenBUGS)

cluster_threshold1 = 3000
cluster_threshold2 = 30000
CLUSTER_START = 1
CLUSTER_STOP = 5
MCMC_NUM = 20000
MCMC_BURN = 2000
MAX_CN = 10
MAX_SCN = 8

gse.id = "GSE42861"

inpath = ""
outpath = ""


#######################
# LOAD + EXECUTE ARGS #
################################################################################################################

# these lines of code are only relevant when running the scripts in batch mode
# they allow to change input parameters before the simulations are executed
# for manually running this script in R, they can be skipped
args = commandArgs(TRUE)
for(i in 1:length(args)){
	eval(parse(text=args[[i]]))
}

set.seed(2021)

source(paste0(inpath,"utility_scripts/utility_functions.R"))
source(paste0(inpath,"utility_scripts/two_arm_model_functions_three_variables.R"))

# depending on your R version either one of these is appropriate and the other won't work: 
#	- "Objects for Deconvolution V1.RData" 
#	- "Objects for Deconvolution V2.RData"
load(paste0(inpath,"deconvolution/Objects for Deconvolution V2.RData"))
source(paste0(inpath,"deconvolution/Deconvolution code.R"))

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 8)
show_col_types = FALSE
gset <- getGEO(GEO = gse.id, filename = NULL)
gset <- gset[[1]]
betas = exprs(gset)
covariates = pData(gset)
covariates = covariates[match(colnames(betas),row.names(covariates)),]
cell_fractions = PredictCellComposition(betas,Array="450K")
cell_fractions[cell_fractions<0] = 0
for(i in 1:nrow(cell_fractions)){
	cell_fractions[i,] = cell_fractions[i,]/sum(cell_fractions[i,])
}

myeloid = rowSums(cell_fractions[,c("Mono","Gran")])
lymphoid = rowSums(cell_fractions[,c("CD4T","CD8T","NK","Bcell")])
cell_fractions = cbind(myeloid,lymphoid)
colnames(cell_fractions) = c("myeloid","lymphoid")

# load annotation data and subset to chromosome
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anchr0 = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anchr0 = as.data.frame(anchr0)
anchr0 = anchr0[row.names(anchr0) %in% row.names(betas),]

gstd_dmcts = read.csv(paste0(inpath,"combined_set.csv"),header=TRUE)
gstd_dmcts = gstd_dmcts$CpG
gstd_dmcts = gstd_dmcts[gstd_dmcts %in% row.names(betas)]
gstd_dmcts = gstd_dmcts[gstd_dmcts %in% row.names(anchr0)]


anno_gstd = anchr0[anchr0$Name %in% gstd_dmcts,]
table(anno_gstd$chr)

TM_gamma_means = list()
TM_gamma_sds = list()
TM_psi_means = list()
TM_psi_sds = list()
TM_tau_means = list()
TM_tau_sds = list()
SCM2_gamma_means = list()
SCM2_gamma_sds = list()
SCM2_psi_means = list()
SCM2_psi_sds = list()
SCM2_tau_means = list()
SCM2_tau_sds = list()
SCM2_ogam_means = list()
SCM2_ogam_sds = list()
SCM2_oogam_means = list()
SCM2_oogam_sds = list()
SCM2_olambda_means = list()
SCM2_olambda_sds = list()
SCM2_oolambda_means = list()
SCM2_oolambda_sds = list()

targetedCpGs = c()
fishedCpGs = c()
gscs_per_cluster = list()

scl.cnt = 0
chr.cnt = 0



chrom = TCHRS

chr.cnt = chr.cnt + 1

anchr = anchr0[anchr0$chr==chrom,]
anchr = anchr[order(anchr$pos),]
betas.chr = betas[row.names(betas) %in% row.names(anchr), ]
betas.chr = betas.chr[match(row.names(anchr),row.names(betas.chr)),]

super_clusters = construct.cluster.structures(
	chr_positions=anchr$pos, cl_threshold1=cluster_threshold1, cl_threshold2=cluster_threshold2,
	max_cn=MAX_CN, max_scn=MAX_SCN, max_split_prop=0.7
)
supercl = super_clusters$sc.v3
cluster_group_length = unlist(lapply(supercl,FUN=length))
supercl = supercl[cluster_group_length > 1]
rm(super_clusters)
gc()

target_clusters = list()
gscs_per_cluster = list()
cnt = 1
cpgs_of_interest = anno_gstd[anno_gstd$chr == chrom,]
scl_position_vector = c()
for( i in 1:length(supercl)){
	scl_position_vector = c(scl_position_vector,unlist(supercl[[i]]))
	if(sum(cpgs_of_interest$pos %in% unlist(supercl[[i]]))>0){
		print(i)
		target_clusters[[cnt]] = supercl[[i]]
		gscs_per_cluster[[cnt]] = sum(cpgs_of_interest$pos %in% unlist(supercl[[i]]))
		cnt = cnt+1
		
	}
}

targetedCpGs = cpgs_of_interest[cpgs_of_interest$pos %in% scl_position_vector,]$Name

print(chrom)

for(rrr in 1:length(target_clusters)){
	scl.cnt = scl.cnt + 1
	print(scl.cnt)
	
	clusters = target_clusters[[rrr]]
	scl_an = anchr[anchr$pos %in% unlist(clusters),]
	
	max_clst_length = max(unlist(lapply(clusters,FUN=length)))
	num_subjects = ncol(betas.chr)
	num_cpgs = c()
	X = array(NA,dim=c(num_subjects,length(clusters),max_clst_length))
	for(clst in 1:length(clusters)){
		TNCPG = length(clusters[[clst]])
		subAN = scl_an[scl_an$pos %in% clusters[[clst]],]
		subBE = t(betas.chr[row.names(betas.chr) %in% row.names(subAN),])
		if(nrow(subAN) == 1){
			subBE = matrix(betas.chr[row.names(betas.chr) %in% row.names(subAN),],ncol=1)
		}
		num_cpgs[clst] = ncol(subBE)
		X[1:nrow(subBE),clst,1:ncol(subBE)] = subBE
	}
	num_ctypes = ncol(cell_fractions)
	
	fishedCpGs = c(fishedCpGs,scl_an[scl_an$pos %in% unlist(clusters),]$Name)
	
	dats = list(
		X = X,
		VX1 = (covariates$"smoking status:ch1"!="never")*1.0,
		VX2 = (covariates$"gender:ch1"=="m")*1.0,
		VX3 = as.numeric(covariates$"age:ch1"),
		N = num_subjects,
		MCPG = num_cpgs,
		DEGCL = (num_cpgs <= 1)*1,
		KCF = num_ctypes,
		CFR = cell_fractions,
		MCL = length(clusters)
	)
	
	inivals = generate.initial.values(clusters=clusters, chains=1, KCF=num_ctypes, N=num_subjects)
	
	print("start fit SCM")
	fit_SCM2_cprior2 <- bugs(
		dats, inits = inivals, 
		parameters.to.save = c(
			"wrapper.mu","wrapper.gam1","tau","psi","omu","oomu",
			"ogam1",
			"oogam1",
			"otau","ootau",
			"olambda1","oolambda1"
		), 
		model.file = SCM2_CP2_3VARS, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
		digits=8
	)
	print("finish fit SCM")
	
	mean_gam_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,9)=="wrapper.g",1],
		nrow=num_ctypes, byrow=TRUE
	)
	sd_gam_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,9)=="wrapper.g",2],
		nrow=num_ctypes, byrow=TRUE
	)
	row.names(mean_gam_SCM2_cp2) = row.names(sd_gam_SCM2_cp2) = colnames(cell_fractions)
	colnames(mean_gam_SCM2_cp2) = colnames(sd_gam_SCM2_cp2) = scl_an[scl_an$pos %in% unlist(clusters),]$Name
	
	mean_ogam_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,4)=="ogam",1],
		nrow=num_ctypes, byrow=TRUE
	)
	sd_ogam_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,4)=="ogam",2],
		nrow=num_ctypes, byrow=TRUE
	)
	row.names(mean_ogam_SCM2_cp2) = row.names(sd_ogam_SCM2_cp2) = colnames(cell_fractions)
	
	mean_oogam_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,5)=="oogam",1]
	sd_oogam_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,5)=="oogam",2]
	names(mean_oogam_SCM2_cp2) = names(sd_oogam_SCM2_cp2) = colnames(cell_fractions)
	
	mean_psi_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="psi",1],
		nrow=num_ctypes, byrow=TRUE
	)
	sd_psi_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="psi",2],
		nrow=num_ctypes, byrow=TRUE
	)
	row.names(mean_psi_SCM2_cp2) = row.names(sd_psi_SCM2_cp2) = colnames(cell_fractions)
	colnames(mean_psi_SCM2_cp2) = colnames(sd_psi_SCM2_cp2) = scl_an[scl_an$pos %in% unlist(clusters),]$Name
	
	mean_olambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,7)=="olambda",1]
	sd_olambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,7)=="olambda",2]
	mean_oolambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,8)=="oolambda",1]
	sd_oolambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,8)=="oolambda",2]
	
	mean_tau_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="tau",1]
	sd_tau_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="tau",2]


	print("start fit TM")
	fit_TM_cprior2 <- bugs(
		dats, inits = inivals, 
		parameters.to.save = c(
			"mu","gam1","tau","psi"
		), 
		model.file = TM_CP2_3VARS, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
		digits=8
	)
	print("finish fit TM")
	
	mean_gam_TM_cp2 = matrix(
		fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="gam",1],
		nrow=num_ctypes, byrow=TRUE
	)
	sd_gam_TM_cp2 = matrix(
		fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="gam",2],
		nrow=num_ctypes, byrow=TRUE
	)
	row.names(mean_gam_TM_cp2) = row.names(sd_gam_TM_cp2) = colnames(cell_fractions)
	colnames(mean_gam_TM_cp2) = colnames(sd_gam_TM_cp2) = scl_an[scl_an$pos %in% unlist(clusters),]$Name
	
	mean_psi_TM_cp2 = matrix(
		fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="psi",1],
		nrow=num_ctypes, byrow=TRUE
	)
	sd_psi_TM_cp2 = matrix(
		fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="psi",2],
		nrow=num_ctypes, byrow=TRUE
	)
	row.names(mean_psi_TM_cp2) = row.names(sd_psi_TM_cp2) = colnames(cell_fractions)
	colnames(mean_psi_TM_cp2) = colnames(sd_psi_TM_cp2) = scl_an[scl_an$pos %in% unlist(clusters),]$Name
	
	mean_tau_TM_cp2 = fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="tau",1]
	sd_tau_TM_cp2 = fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="tau",2]
	

	TM_gamma_means[[rrr]] = mean_gam_TM_cp2
	TM_gamma_sds[[rrr]] = sd_gam_TM_cp2
	TM_psi_means[[rrr]] = mean_psi_TM_cp2
	TM_psi_sds[[rrr]] = sd_psi_TM_cp2
	TM_tau_means[[rrr]] = mean_tau_TM_cp2
	TM_tau_sds[[rrr]] = sd_tau_TM_cp2
	SCM2_gamma_means[[rrr]] = mean_gam_SCM2_cp2
	SCM2_gamma_sds[[rrr]] = sd_gam_SCM2_cp2
	SCM2_psi_means[[rrr]] = mean_psi_SCM2_cp2
	SCM2_psi_sds[[rrr]] = sd_psi_SCM2_cp2
	SCM2_tau_means[[rrr]] = mean_tau_SCM2_cp2
	SCM2_tau_sds[[rrr]] = sd_tau_SCM2_cp2
	SCM2_ogam_means[[rrr]] = mean_ogam_SCM2_cp2
	SCM2_ogam_sds[[rrr]] = sd_ogam_SCM2_cp2
	SCM2_oogam_means[[rrr]] = mean_oogam_SCM2_cp2
	SCM2_oogam_sds[[rrr]] = sd_oogam_SCM2_cp2
	SCM2_olambda_means[[rrr]] = mean_olambda_SCM2_cp2
	SCM2_olambda_sds[[rrr]] = sd_olambda_SCM2_cp2
	SCM2_oolambda_means[[rrr]] = mean_oolambda_SCM2_cp2
	SCM2_oolambda_sds[[rrr]] = sd_oolambda_SCM2_cp2
	
	out_string = paste0(
		outpath,
		"run.",gse.id,"_cth1.",cluster_threshold1,"_cth2.",cluster_threshold2,"_",chrom,"_mat.RDATA"
	)
	save(
		TM_gamma_means, TM_gamma_sds, TM_psi_means, TM_psi_sds, TM_tau_means, TM_tau_sds,
		SCM2_gamma_means, SCM2_gamma_sds, SCM2_psi_means, SCM2_psi_sds, SCM2_tau_means, SCM2_tau_sds,
		SCM2_ogam_means, SCM2_ogam_sds, SCM2_oogam_means, SCM2_oogam_sds,
		SCM2_olambda_means,SCM2_olambda_sds,SCM2_oolambda_means,SCM2_oolambda_sds,
		fishedCpGs,targetedCpGs,gscs_per_cluster,target_clusters,
		file=out_string
	)
	cat("Saved "); print(out_string)
	
	logt = data.frame(chromosomes=chr.cnt,sclusts=scl.cnt,current=chrom)
	write.csv(logt,file=paste0(outpath,"log.",chrom,".",gse.id,".txt"))
}

