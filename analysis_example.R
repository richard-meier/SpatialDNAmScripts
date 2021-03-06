
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

cluster_threshold1 = 3000
cluster_threshold2 = 30000
MCMC_NUM = 15000
MCMC_BURN = 1000
MAX_CN = 10
MAX_SCN = 8


# specify the project directory:
inpath = ""


source(paste0(inpath,"utility_scripts/utility_functions.R"))
source(paste0(inpath,"utility_scripts/two_arm_model_functions_analysis.R"))


# load dataset (beta values, covariates, cell fractions)
load(paste0(inpath,"data/example_data.RData"))


# load annotation data and subset chromosome 1 loci to those shared between annotation and dataset
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anchr = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anchr = as.data.frame(anchr)
anchr = anchr[row.names(anchr) %in% row.names(betas),]
anchr = anchr[anchr$chr=="chr1",]
anchr = anchr[order(anchr$pos),]
betas.chr = betas[row.names(betas) %in% row.names(anchr), ]
betas.chr = betas.chr[match(row.names(anchr),row.names(betas.chr)),]


# perform cluster generating algorithm
super_clusters = construct.cluster.structures(
	chr_positions=anchr$pos, cl_threshold1=cluster_threshold1, cl_threshold2=cluster_threshold2,
	max_cn=MAX_CN, max_scn=MAX_SCN, max_split_prop=0.7
)
supercl = super_clusters$sc.v3
cluster_group_length = unlist(lapply(supercl,FUN=length))
supercl = supercl[cluster_group_length > 1]
rm(super_clusters); gc()

length(supercl)
length(unlist(supercl[[1]]))
length(unlist(supercl[[2]]))
length(unlist(supercl[[3]]))
length(unlist(supercl[[4]]))
length(unlist(supercl[[5]]))



# create lists to store output
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

# run analysis of the whole dataset. NOTE: this will take a while!
nCpGs = 0
t_start = Sys.time()
for(cid in 1:length(supercl)){
	clusters = supercl[[cid]]
	
	cat("Processing cluster ",cid," containing ",length(unlist(clusters))," CpGs ...\n",sep="")
	
	
	# restructure data for the model fit
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
		GRP = covariates$case,
		N = num_subjects,
		MCPG = num_cpgs,
		DEGCL = (num_cpgs <= 1)*1,
		KCF = ncol(cell_fractions),
		CFR = cell_fractions,
		MCL = length(clusters)
	)
	
	
	# generate initial values for the MCMC sampler
	inivals = generate.initial.values(clusters=clusters, chains=1, KCF=ncol(cell_fractions), N=num_subjects)


	# fit the scm2 model and store its output
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


	# fit the TM model and store its output
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
	

	TM_gamma_means[[cid]] = mean_gam_TM_cp2
	TM_gamma_sds[[cid]] = sd_gam_TM_cp2
	TM_psi_means[[cid]] = mean_psi_TM_cp2
	TM_psi_sds[[cid]] = sd_psi_TM_cp2
	TM_tau_means[[cid]] = mean_tau_TM_cp2
	TM_tau_sds[[cid]] = sd_tau_TM_cp2
	SCM2_gamma_means[[cid]] = mean_gam_SCM2_cp2
	SCM2_gamma_sds[[cid]] = sd_gam_SCM2_cp2
	SCM2_psi_means[[cid]] = mean_psi_SCM2_cp2
	SCM2_psi_sds[[cid]] = sd_psi_SCM2_cp2
	SCM2_tau_means[[cid]] = mean_tau_SCM2_cp2
	SCM2_tau_sds[[cid]] = sd_tau_SCM2_cp2
	SCM2_ogam_means[[cid]] = mean_ogam_SCM2_cp2
	SCM2_ogam_sds[[cid]] = sd_ogam_SCM2_cp2
	SCM2_oogam_means[[cid]] = mean_oogam_SCM2_cp2
	SCM2_oogam_sds[[cid]] = sd_oogam_SCM2_cp2
	SCM2_olambda_means[[cid]] = mean_olambda_SCM2_cp2
	SCM2_olambda_sds[[cid]] = sd_olambda_SCM2_cp2
	SCM2_oolambda_means[[cid]] = mean_oolambda_SCM2_cp2
	SCM2_oolambda_sds[[cid]] = sd_oolambda_SCM2_cp2
	
	nCpGs = nCpGs + length(unlist(clusters))
	cat("Finished ",cid, " super clusters and ", nCpGs," CpGs after a ",sep="")
	print(Sys.time()-t_start)
	cat("\n")
}


# assemble model estimates into matrices
gam_means_SCM2 = matrix(NA,nrow=6,ncol=0)
gam_sds_SCM2 = matrix(NA,nrow=6,ncol=0)
tau_means_mat_SCM2 = matrix(NA,nrow=6,ncol=0)
gam_means_TM = matrix(NA,nrow=6,ncol=0)
gam_sds_TM = matrix(NA,nrow=6,ncol=0)
tau_means_mat_TM = matrix(NA,nrow=6,ncol=0)
for(cid in 1:length(supercl)){
	gam_means_SCM2 = cbind(gam_means_SCM2, SCM2_gamma_means[[cid]])
	gam_sds_SCM2 = cbind(gam_sds_SCM2, SCM2_gamma_sds[[cid]])
	tau_tmp = SCM2_psi_means[[cid]]*0 + SCM2_tau_means[[cid]]
	tau_means_mat_SCM2 = cbind(tau_means_mat_SCM2, tau_tmp)
	
	gam_means_TM = cbind(gam_means_TM, TM_gamma_means[[cid]])
	gam_sds_TM = cbind(gam_sds_TM, TM_gamma_sds[[cid]])
	tau_tmp = TM_psi_means[[cid]]*0 + TM_tau_means[[cid]]
	tau_means_mat_TM = cbind(tau_means_mat_TM, tau_tmp)
}


# calculate kappa factors
kf_mat_SCM2 = sqrt(1/tau_means_mat_SCM2)
kf_mat_SCM2 = prfun(x=kf_mat_SCM2,0.6950,8.0301,-25.9266,minval=0.0,maxval=0.1128)
kf_mat_TM = sqrt(1/tau_means_mat_TM)
kf_mat_TM = kf_mat_TM = prfun(x=kf_mat_TM,0.8525,12.9091,-52.1179,minval=0.0,maxval=0.1143)


# calculate tests statistics (pseudo-p-values)
SCM2_bpv_mat = gam_means_SCM2 / (gam_sds_SCM2*kf_mat_SCM2)
for(i in 1:6){
	SCM2_bpv_mat[i,] = (pnorm(q=-abs(SCM2_bpv_mat[i,])))*2
}
TM_bpv_mat = gam_means_TM / (gam_sds_TM*kf_mat_TM)
for(i in 1:6){
	TM_bpv_mat[i,] = (pnorm(q=-abs(TM_bpv_mat[i,])))*2
}


# create SCM2 result table
i=1
results_scm2 = data.frame(
	loci=colnames(SCM2_bpv_mat), cell.type=row.names(SCM2_bpv_mat)[i],
	delta.mean=gam_means_SCM2[i,], delta.sd=gam_sds_SCM2[i,], 
	kappa=kf_mat_SCM2[i,], test.stat=SCM2_bpv_mat[i,], 
	reject.h0=SCM2_bpv_mat[i,]<0.05
)
for(i in 2:6){
	res_tmp = data.frame(
		loci=colnames(SCM2_bpv_mat), cell.type=row.names(SCM2_bpv_mat)[i],
		delta.mean=gam_means_SCM2[i,], delta.sd=gam_sds_SCM2[i,], 
		kappa=kf_mat_SCM2[i,], test.stat=SCM2_bpv_mat[i,], 
		reject.h0=SCM2_bpv_mat[i,]<0.05
	)
	results_scm2 = rbind( results_scm2, res_tmp )
}


# create TM result table
i=1
results_tm = data.frame(
	loci=colnames(TM_bpv_mat), cell.type=row.names(TM_bpv_mat)[i],
	delta.mean=gam_means_TM[i,], delta.sd=gam_sds_TM[i,], 
	kappa=kf_mat_TM[i,], test.stat=TM_bpv_mat[i,], 
	reject.h0=TM_bpv_mat[i,]<0.05
)
for(i in 2:6){
	res_tmp = data.frame(
		loci=colnames(TM_bpv_mat), cell.type=row.names(TM_bpv_mat)[i],
		delta.mean=gam_means_TM[i,], delta.sd=gam_sds_TM[i,], 
		kappa=kf_mat_TM[i,], test.stat=TM_bpv_mat[i,], 
		reject.h0=TM_bpv_mat[i,]<0.05
	)
	results_tm = rbind( results_tm, res_tmp )
}

# print first 5 entries in result tables
head(results_scm2, n=5)
head(results_tm, n=5)
