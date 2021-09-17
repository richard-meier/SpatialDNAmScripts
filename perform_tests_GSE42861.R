
options(stringsAsFactors=FALSE)

library(gtools)
library(MASS)
library(EpiDISH)
library(TCA) # requires version 1.2 !
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GEOquery)

gse.id = "GSE42861"
CLIMIT = 136

# specify the project directory:
inpath = ""

# specify output directory where files will be saved:
outpath = ""

source(paste0(inpath,"utility_scripts/utility_functions.R"))


# depending on your R version either one of these is approprite and the other won't work: 
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


# load target list of CpGs
gstd_dmct_tab = read.csv(paste0(inpath,"data/combined_set.csv"),header=TRUE)
gstd_dmcts = gstd_dmct_tab$CpG


# load SCM2 and TM2 model fits
tdir = paste0(inpath,"data_analysis/run.GSE42861_cth1.3000_cth2.30000_")
gam_means_SCM2 = matrix(NA,nrow=2,ncol=0)
gam_sds_SCM2 = matrix(NA,nrow=2,ncol=0)
tau_means_mat_SCM2 = matrix(NA,nrow=2,ncol=0)
gam_means_TM = matrix(NA,nrow=2,ncol=0)
gam_sds_TM = matrix(NA,nrow=2,ncol=0)
tau_means_mat_TM = matrix(NA,nrow=2,ncol=0)
chroms = c(paste0("chr",1:22),"chrX")
for(chr in chroms){
	load(paste0(tdir,chr,"_mat.RDATA"))
	for( r in 1:min(length(SCM2_gamma_means),CLIMIT) ){
		gam_means_SCM2 = cbind(gam_means_SCM2, SCM2_gamma_means[[r]])
		gam_sds_SCM2 = cbind(gam_sds_SCM2, SCM2_gamma_sds[[r]])
		tau_tmp = SCM2_psi_means[[r]]*0 + SCM2_tau_means[[r]]
		tau_means_mat_SCM2 = cbind(tau_means_mat_SCM2, tau_tmp)
		
		gam_means_TM = cbind(gam_means_TM, TM_gamma_means[[r]])
		gam_sds_TM = cbind(gam_sds_TM, TM_gamma_sds[[r]])
		tau_tmp = TM_psi_means[[r]]*0 + TM_tau_means[[r]]
		tau_means_mat_TM = cbind(tau_means_mat_TM, tau_tmp)
	}
}
processed_targets = colnames(gam_means_SCM2) %in% gstd_dmcts


# calculate kappa factors
kf_mat_SCM2 = sqrt(1/tau_means_mat_SCM2)
kf_mat_SCM2 = prfun(x=kf_mat_SCM2,0.6950,8.0301,-25.9266,minval=0.0,maxval=0.1128)
kf_mat_TM = sqrt(1/tau_means_mat_TM)
kf_mat_TM = prfun(x=kf_mat_TM,0.8525,12.9091,-52.1179,minval=0.0,maxval=0.1143)


# calculate test statistics (pseudo-p-values)
SCM2_bpv_mat = gam_means_SCM2 / (gam_sds_SCM2*kf_mat_SCM2)
for(i in 1:2){
	SCM2_bpv_mat[i,] = (pnorm(q=-abs(SCM2_bpv_mat[i,])))*2
}
TM_bpv_mat = gam_means_TM / (gam_sds_TM*kf_mat_TM)
for(i in 1:2){
	TM_bpv_mat[i,] = (pnorm(q=-abs(TM_bpv_mat[i,])))*2
}


# prepare Bayesian model results
tests_SCM2_myeloid = data.frame(
	Estimate = gam_means_SCM2["myeloid",],
	p = SCM2_bpv_mat["myeloid",],
	row.names = colnames(gam_means_SCM2)
)
tests_SCM2_lymphoid = data.frame(
	Estimate = gam_means_SCM2["lymphoid",],
	p = SCM2_bpv_mat["lymphoid",],
	row.names = colnames(gam_means_SCM2)
)
tests_TM_myeloid = data.frame(
	Estimate = gam_means_TM["myeloid",],
	p = TM_bpv_mat["myeloid",],
	row.names = colnames(gam_means_TM)
)
tests_TM_lymphoid = data.frame(
	Estimate = gam_means_TM["lymphoid",],
	p = TM_bpv_mat["lymphoid",],
	row.names = colnames(gam_means_TM)
)


# load annotation data and subset to chromosome
betas.chr = betas[row.names(betas) %in% row.names(tests_SCM2_myeloid), ]
betas.chr = betas.chr[match(row.names(tests_SCM2_myeloid),row.names(betas.chr)),]


# process data with CellDMC
smoker = (covariates$"smoking status:ch1"!="never")*1.0
covs = cbind(
	(covariates$"gender:ch1"=="m")*1.0,
	as.numeric(covariates$"age:ch1")
)
colnames(covs) = c("gender","age")
cstart = Sys.time()
fit_cdmc <- CellDMC(beta.m=betas.chr, pheno.v=smoker, frac.m=cell_fractions,cov.mod=covs)
tests_cdmc_myeloid = fit_cdmc$coe[["myeloid"]]
tests_cdmc_lymphoid = fit_cdmc$coe[["lymphoid"]]
print(Sys.time()-cstart)


# process data with TCA
C1 = cbind(smoker,covs)
row.names(C1) = row.names(cell_fractions)
cstart = Sys.time()
fit_tca = tca(
	X = betas.chr, W = cell_fractions, C1 = C1, C2 = NULL,
	parallel = TRUE, num_cores = 1, verbose = FALSE, 
	log_file = paste0(inpath,"tca_log.txt")
)
tests_tca_myeloid = data.frame(
	Estimate = fit_tca$gammas_hat[,"myeloid.smoker"],
	p = fit_tca$gammas_hat_pvals[,"myeloid.smoker"]
)
tests_tca_lymphoid = data.frame(
	Estimate = fit_tca$gammas_hat[,"lymphoid.smoker"],
	p = fit_tca$gammas_hat_pvals[,"lymphoid.smoker"]
)
print(Sys.time()-cstart)


# compile result tables
full_tab_myeloid = data.frame(
	loci = row.names(tests_SCM2_myeloid),
	delta.scm2 = tests_SCM2_myeloid$Estimate,
	p.scm2 = tests_SCM2_myeloid$p,
	delta.tm = tests_TM_myeloid$Estimate,
	p.tm = tests_TM_myeloid$p,
	delta.cdmc = tests_cdmc_myeloid$Estimate,
	p.cdmc = tests_cdmc_myeloid$p,
	delta.tca = tests_tca_myeloid$Estimate,
	p.tca = tests_tca_myeloid$p
)
full_tab_lymphoid = data.frame(
	loci = row.names(tests_SCM2_lymphoid),
	delta.scm2 = tests_SCM2_lymphoid$Estimate,
	p.scm2 = tests_SCM2_lymphoid$p,
	delta.tm = tests_TM_lymphoid$Estimate,
	p.tm = tests_TM_lymphoid$p,
	delta.cdmc = tests_cdmc_lymphoid$Estimate,
	p.cdmc = tests_cdmc_lymphoid$p,
	delta.tca = tests_tca_lymphoid$Estimate,
	p.tca = tests_tca_lymphoid$p
)
write.csv(full_tab_myeloid,file=paste0(outpath,"tests_myeloid.csv"),row.names=FALSE)
write.csv(full_tab_lymphoid,file=paste0(outpath,"tests_lymphoid.csv"),row.names=FALSE)

