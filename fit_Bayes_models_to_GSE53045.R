
install.packages("gtools")
install.packages("R2OpenBUGS")
install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

options(stringsAsFactors=FALSE)

library(GEOquery)
library(gtools)
library(R2OpenBUGS)
library(MASS)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

chrom = "chr1"
cluster_threshold1 = 3000
cluster_threshold2 = 30000
CLUSTER_START = 1
CLUSTER_STOP = 5
MCMC_NUM = 20000
MCMC_BURN = 2000
MAX_CN = 10
MAX_SCN = 8

# specify the project directory:
inpath = ""

# specify output directory where files will be saved:
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


################
# PROGRAM CODE #
################################################################################################################

gse.id = "GSE53045"
set.seed(2021)

source(paste0(inpath,"utility_scripts/utility_functions.R"))
source(paste0(inpath,"utility_scripts/two_arm_model_functions_analysis.R"))

# depending on your R version either one of these is approprite and the other won't work: 
#	- "Objects for Deconvolution V1.RData" 
#	- "Objects for Deconvolution V2.RData"
load(paste0(inpath,"deconvolution/Objects for Deconvolution V2.RData"))
source(paste0(inpath,"deconvolution/Deconvolution code.R"))

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

# load annotation data and subset to chromosome
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anchr = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anchr = as.data.frame(anchr)
anchr = anchr[row.names(anchr) %in% row.names(betas),]
anchr = anchr[anchr$chr==chrom,]
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

nCpGs = 0
rrr = 1
t_start = Sys.time()
CLUSTER_STOP = min(CLUSTER_STOP,length(supercl))
for(rrr in CLUSTER_START:CLUSTER_STOP){
	clusters = supercl[[rrr]]
	
	cat("Processing cluster ",rrr," containing ",length(unlist(clusters))," CpGs ...\n",sep="")
	
	max_clst_length = max(unlist(lapply(clusters,FUN=length)))
	num_subjects = ncol(betas.chr)
	num_cpgs = c()
	X = array(NA,dim=c(num_subjects,length(clusters),max_clst_length))
	for(clst in 1:length(clusters)){
		TNCPG = length(clusters[[clst]])
		subAN = anchr[anchr$pos %in% clusters[[clst]],]
		subBE = t(betas.chr[row.names(betas.chr) %in% row.names(subAN),])
		if(nrow(subAN) == 1){
			subBE = matrix(betas.chr[row.names(betas.chr) %in% row.names(subAN),],ncol=1)
		}
		num_cpgs[clst] = ncol(subBE)
		X[1:nrow(subBE),clst,1:ncol(subBE)] = subBE
	}
	
	dats = list(
		X = X,
		GRP = (covariates$"disease state:ch1"=="Smoker")*1.0,
		N = num_subjects,
		MCPG = num_cpgs,
		DEGCL = (num_cpgs <= 1)*1,
		KCF = ncol(cell_fractions),
		CFR = cell_fractions,
		MCL = length(clusters)
	)
	
	
	inivals = generate.initial.values(clusters=clusters, chains=1, KCF=6, N=num_subjects)

	fit_SCM2_cprior2 <- bugs(
		dats, inits = inivals, 
		parameters.to.save = c(
			"wrapper.mu","wrapper.gam","tau","psi","omu","oomu","ogam","oogam",
			"otau","ootau","olambda","oolambda"
		), 
		model.file = SCM2_CP2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
		digits=8 #, debug=TRUE
	)
	
	mean_gam_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,9)=="wrapper.g",1],
		nrow=6, byrow=TRUE
	)
	sd_gam_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,9)=="wrapper.g",2],
		nrow=6, byrow=TRUE
	)
	row.names(mean_gam_SCM2_cp2) = row.names(sd_gam_SCM2_cp2) = colnames(cell_fractions)
	colnames(mean_gam_SCM2_cp2) = colnames(sd_gam_SCM2_cp2) = anchr[anchr$pos %in% unlist(clusters),]$Name
	
	mean_oogam_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,5)=="oogam",1]
	sd_oogam_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,5)=="oogam",2]
	names(mean_oogam_SCM2_cp2) = names(sd_oogam_SCM2_cp2) = colnames(cell_fractions)
	
	mean_ogam_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,4)=="ogam",1],
		nrow=6, byrow=TRUE
	)
	sd_ogam_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,4)=="ogam",2],
		nrow=6, byrow=TRUE
	)
	row.names(mean_ogam_SCM2_cp2) = row.names(sd_ogam_SCM2_cp2) = colnames(cell_fractions)
	
	mean_oogam_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,5)=="oogam",1]
	sd_oogam_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,5)=="oogam",2]
	names(mean_oogam_SCM2_cp2) = names(sd_oogam_SCM2_cp2) = colnames(cell_fractions)
	
	mean_psi_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="psi",1],
		nrow=6, byrow=TRUE
	)
	sd_psi_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="psi",2],
		nrow=6, byrow=TRUE
	)
	row.names(mean_psi_SCM2_cp2) = row.names(sd_psi_SCM2_cp2) = colnames(cell_fractions)
	colnames(mean_psi_SCM2_cp2) = colnames(sd_psi_SCM2_cp2) = anchr[anchr$pos %in% unlist(clusters),]$Name
	
	mean_olambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,7)=="olambda",1]
	sd_olambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,7)=="olambda",2]
	mean_oolambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,8)=="oolambda",1]
	sd_oolambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,8)=="oolambda",2]
	
	mean_tau_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="tau",1]
	sd_tau_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="tau",2]


	fit_TM_cprior2 <- bugs(
		dats, inits = inivals, 
		parameters.to.save = c(
			"mu","gam","tau","psi"
		), 
		model.file = TM_CP2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
		digits=8 #, debug=TRUE
	)
	
	mean_gam_TM_cp2 = matrix(
		fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,4)=="gam[",1],
		nrow=6, byrow=TRUE
	)
	sd_gam_TM_cp2 = matrix(
		fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,4)=="gam[",2],
		nrow=6, byrow=TRUE
	)
	row.names(mean_gam_TM_cp2) = row.names(sd_gam_TM_cp2) = colnames(cell_fractions)
	colnames(mean_gam_TM_cp2) = colnames(sd_gam_TM_cp2) = anchr[anchr$pos %in% unlist(clusters),]$Name
	
	mean_psi_TM_cp2 = matrix(
		fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="psi",1],
		nrow=6, byrow=TRUE
	)
	sd_psi_TM_cp2 = matrix(
		fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="psi",2],
		nrow=6, byrow=TRUE
	)
	row.names(mean_psi_TM_cp2) = row.names(sd_psi_TM_cp2) = colnames(cell_fractions)
	colnames(mean_psi_TM_cp2) = colnames(sd_psi_TM_cp2) = anchr[anchr$pos %in% unlist(clusters),]$Name
	
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
		"run.",gse.id,"_cth1.",cluster_threshold1,"_cth2.",cluster_threshold2,"_",chrom, "_start.",
		CLUSTER_START,"_stop.",CLUSTER_STOP,"_mat.RDATA"
	)
	save(
		TM_gamma_means, TM_gamma_sds, TM_psi_means, TM_psi_sds, TM_tau_means, TM_tau_sds,
		SCM2_gamma_means, SCM2_gamma_sds, SCM2_psi_means, SCM2_psi_sds, SCM2_tau_means, SCM2_tau_sds,
		SCM2_ogam_means, SCM2_ogam_sds, SCM2_oogam_means, SCM2_oogam_sds,
		SCM2_olambda_means,SCM2_olambda_sds,SCM2_oolambda_means,SCM2_oolambda_sds,
		file=out_string
	)
	
	nCpGs = nCpGs + length(unlist(clusters))
	cat("Finished ",rrr, " super clusters and ", nCpGs," CpGs after a ",sep="")
	print(Sys.time()-t_start)
}
