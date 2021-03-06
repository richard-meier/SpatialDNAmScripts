
install.packages("gtools")
install.packages("TCA")
install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("EpiDISH")

options(stringsAsFactors=FALSE)

library(gtools)
library(MASS)
library(EpiDISH)
library(TCA) # requires version 1.2 !
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


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
betas.chr = na.omit(betas.chr)


# process data with CellDMC
cstart = Sys.time()
smoker = (covariates$"disease state:ch1"=="Smoker")*1.0
fit_cdmc <- CellDMC(beta.m=betas.chr, pheno.v=smoker, frac.m=cell_fractions)
tests_cdmc = fit_cdmc$coe[["CD8T"]]
print(Sys.time()-cstart)


# process data with TCA
cstart = Sys.time()
C1 = matrix(smoker,ncol=1)
colnames(C1)="smoker"
row.names(C1) = row.names(cell_fractions)
fit_tca = tca(
	X = betas.chr, W = cell_fractions, C1 = C1, C2 = NULL,
	parallel = TRUE, num_cores = 8, verbose = FALSE, 
	log_file = paste0(inpath,"tca_log.txt")
)
tests_tca = fit_tca$gammas_hat
pvs_tca = fit_tca$gammas_hat_pvals
print(Sys.time()-cstart)


# load the 26 model fit batches in which chromosome 1 was processed with the Bayes Models
# these batches on the "fit_Bayes_models_to_GSE53045.R" script were run as follows:
# from    1 to  100
# from  101 to  200
#  ... ... ... ...
# from 2401 to 2500
# from 2501 to 2564
gam_means_SCM2 = matrix(NA,nrow=6,ncol=0)
gam_sds_SCM2 = matrix(NA,nrow=6,ncol=0)
tau_means_mat_SCM2 = matrix(NA,nrow=6,ncol=0)
gam_means_TM = matrix(NA,nrow=6,ncol=0)
gam_sds_TM = matrix(NA,nrow=6,ncol=0)
tau_means_mat_TM = matrix(NA,nrow=6,ncol=0)
batches_start = 1+(0:25)*100
batches_stop = (1:26)*100
batches_stop[26] = 2564
for(b in 1:26){
	load(paste0(
		"S:/Biostats/BIO-STAT/Koestler Devin/GRAs/Richard Meier/Project3/real_data_analysis/",
		"run.GSE53045_cth1.3000_cth2.30000_chr1_start.",batches_start[b],
		"_stop.",batches_stop[b],
		"_mat.RDATA"
	))
	for( r in batches_start[b]:batches_stop[b] ){
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


# this needs to be updated according to where the output was stored after running "processGSE147430.R"
load("...\\GSE147430_testing.RDATA")


# identify common CpGs shared between datasets
common_cpgs = intersect(row.names(tests_cdmc), colnames(gam_means_SCM2))
common_cpgs = intersect(common_cpgs, colnames(gam_means_TM))
common_cpgs = intersect(common_cpgs, row.names(tests_tca))
common_cpgs = intersect(common_cpgs, row.names(isolate_tests))


# subset Bayesian model results and isolate results to common CpGs
gset_SCM2 = gam_means_SCM2[,colnames(gam_means_SCM2) %in% common_cpgs]
sdset_SCM2 = gam_sds_SCM2[,colnames(gam_sds_SCM2) %in% common_cpgs]
bpvset_SCM2 = SCM2_bpv_mat[,colnames(SCM2_bpv_mat) %in% common_cpgs]
gset_TM = gam_means_TM[,colnames(gam_means_TM) %in% common_cpgs]
sdset_TM = gam_sds_TM[,colnames(gam_sds_TM) %in% common_cpgs]
bpvset_TM = TM_bpv_mat[,colnames(TM_bpv_mat) %in% common_cpgs]
iso = isolate_tests[row.names(isolate_tests) %in% common_cpgs,]
iso = iso[match(colnames(gset_SCM2),row.names(iso)),]


# bring CellDMC output into the same form as Bayesian models
pvs_cdmc = bpvset_SCM2 * NA
gset_cdmc = bpvset_SCM2 * NA
for(i in 1:nrow(bpvset_SCM2)){
	ctype = row.names(bpvset_SCM2)[i]
	tmp_tab = fit_cdmc$coe[[ctype]]
	tmp_tab = tmp_tab[row.names(tmp_tab) %in% common_cpgs,]
	tmp_tab = tmp_tab[match(colnames(bpvset_SCM2),row.names(tmp_tab)),]
	pvs_cdmc[i,] = tmp_tab$p 
	gset_cdmc[i,] = tmp_tab$Estimate
}

# bring TCA output into the same form as Bayesian models
pvs_tca = t(fit_tca$gammas_hat_pvals)
pvs_tca = pvs_tca[,colnames(pvs_tca) %in% colnames(bpvset_SCM2)]
pvs_tca = pvs_tca[,match(colnames(bpvset_SCM2),colnames(pvs_tca))]
gset_tca = t(fit_tca$gammas_hat)
gset_tca = gset_tca[,colnames(gset_tca) %in% colnames(bpvset_SCM2)]
gset_tca = gset_tca[,match(colnames(bpvset_SCM2),colnames(gset_tca))]
row.names(pvs_tca) = row.names(gset_tca) = row.names(gset_cdmc)


# create full result table
iso_sub = iso
gamsub_scm2 = gset_SCM2[,colnames(gset_SCM2) %in% row.names(iso_sub)]
gamsub_tm = gset_TM[,colnames(gset_TM) %in% row.names(iso_sub)]
gamsub_cdmc = gset_cdmc[,colnames(gset_cdmc) %in% row.names(iso_sub)]
gamsub_tca = gset_tca[,colnames(gset_tca) %in% row.names(iso_sub)]
pvsub_scm2 = SCM2_bpv_mat[,colnames(SCM2_bpv_mat) %in% row.names(iso_sub)]
pvsub_tm = TM_bpv_mat[,colnames(TM_bpv_mat) %in% row.names(iso_sub)]
pvsub_cdmc = pvs_cdmc[,colnames(pvs_cdmc) %in% row.names(iso_sub)]
pvsub_tca = pvs_tca[,colnames(pvs_tca) %in% row.names(iso_sub)]
full_tab = data.frame(
	loci=row.names(iso_sub),
	delta.iso=iso_sub$Estimate,
	p.iso=iso_sub$Pr,
	delta.scm2=gamsub_scm2[2,],
	pseudo.p.scm2=pvsub_scm2[2,],
	delta.tm=gamsub_tm[2,],
	pseudo.p.tm=pvsub_tm[2,],
	delta.cdmc=gamsub_cdmc[2,],
	pseudo.p.cdmc=pvsub_cdmc[2,],
	delta.tca=gamsub_tca[2,],
	pseudo.p.tca=pvsub_tca[2,]
)
row.names(full_tab) = NULL
write.csv(sub_tab,file=paste0(outpath,"tests_full.csv"),row.names=FALSE)
pdf(file=paste0(outpath,"point_estimates_full.pdf"),width=11,height=3)
par(mfrow=c(1,4))
plot(x=iso_sub$Estimate,y=gamsub_scm2[2,],ylim=c(-1,1),xlim=c(-0.15,0.15),col=rgb(0,0,0,0.34),main="SCM2",xlab="isolate delta",ylab="predicted delta"); lines(c(-9,9),c(-9,9),col="red")
plot(x=iso_sub$Estimate,y=gamsub_tm[2,],ylim=c(-2.5,2.5),xlim=c(-0.15,0.15),col=rgb(0,0,0,0.34),main="TM",xlab="isolate delta",ylab="predicted delta"); lines(c(-9,9),c(-9,9),col="red")
plot(x=iso_sub$Estimate,y=gamsub_tca[2,],ylim=c(-2.5,2.5),xlim=c(-0.15,0.15),col=rgb(0,0,0,0.34),main="TCA",xlab="isolate delta",ylab="predicted delta"); lines(c(-9,9),c(-9,9),col="red")
plot(x=iso_sub$Estimate,y=gamsub_cdmc[2,],ylim=c(-1,1),xlim=c(-0.15,0.15),col=rgb(0,0,0,0.34),main="CellDMC",xlab="isolate delta",ylab="predicted delta"); lines(c(-9,9),c(-9,9),col="red")
dev.off()


# create subset result table
iso_sub = iso[abs(iso$Estimate) > 0.05,]
iso_sub = iso_sub[iso_sub$Pr < 0.01,]
gamsub_scm2 = gset_SCM2[,colnames(gset_SCM2) %in% row.names(iso_sub)]
gamsub_tm = gset_TM[,colnames(gset_TM) %in% row.names(iso_sub)]
gamsub_cdmc = gset_cdmc[,colnames(gset_cdmc) %in% row.names(iso_sub)]
gamsub_tca = gset_tca[,colnames(gset_tca) %in% row.names(iso_sub)]
pvsub_scm2 = SCM2_bpv_mat[,colnames(SCM2_bpv_mat) %in% row.names(iso_sub)]
pvsub_tm = TM_bpv_mat[,colnames(TM_bpv_mat) %in% row.names(iso_sub)]
pvsub_cdmc = pvs_cdmc[,colnames(pvs_cdmc) %in% row.names(iso_sub)]
pvsub_tca = pvs_tca[,colnames(pvs_tca) %in% row.names(iso_sub)]
sub_tab = data.frame(
	loci=row.names(iso_sub),
	delta.iso=iso_sub$Estimate,
	p.iso=iso_sub$Pr,
	delta.scm2=gamsub_scm2[2,],
	pseudo.p.scm2=pvsub_scm2[2,],
	delta.tm=gamsub_tm[2,],
	pseudo.p.tm=pvsub_tm[2,],
	delta.cdmc=gamsub_cdmc[2,],
	pseudo.p.cdmc=pvsub_cdmc[2,],
	delta.tca=gamsub_tca[2,],
	pseudo.p.tca=pvsub_tca[2,]
)
row.names(sub_tab) = NULL
write.csv(sub_tab,file=paste0(outpath,"tests_subset.csv"),row.names=FALSE)
pdf(file=paste0(outpath,"point_estimates_subset.pdf"),width=11,height=3)
par(mfrow=c(1,4))
plot(x=iso_sub$Estimate,y=gamsub_scm2[2,],ylim=c(-1,1),xlim=c(-0.15,0.15),col=rgb(0,0,0,0.34),main="SCM2",xlab="isolate delta",ylab="predicted delta"); lines(c(-9,9),c(-9,9),col="red")
plot(x=iso_sub$Estimate,y=gamsub_tm[2,],ylim=c(-1,1),xlim=c(-0.15,0.15),col=rgb(0,0,0,0.34),main="TM",xlab="isolate delta",ylab="predicted delta"); lines(c(-9,9),c(-9,9),col="red")
plot(x=iso_sub$Estimate,y=gamsub_tca[2,],ylim=c(-1,1),xlim=c(-0.15,0.15),col=rgb(0,0,0,0.34),main="TCA",xlab="isolate delta",ylab="predicted delta"); lines(c(-9,9),c(-9,9),col="red")
plot(x=iso_sub$Estimate,y=gamsub_cdmc[2,],ylim=c(-1,1),xlim=c(-0.15,0.15),col=rgb(0,0,0,0.34),main="CellDMC",xlab="isolate delta",ylab="predicted delta"); lines(c(-9,9),c(-9,9),col="red")
dev.off()


# check recovery metrics

orientation_consistent_scm2 = (((gamsub_scm2[2,]<0) & (iso_sub$Estimate<0)) | ((gamsub_scm2[2,]>0) & (iso_sub$Estimate>0)))
significant_consistent_scm2 = pvsub_scm2[2,]<0.05
round(c(
sqrt(mean((iso_sub$Estimate - gamsub_scm2[2,])^2)),
mean(significant_consistent_scm2) * 100,
mean(orientation_consistent_scm2) * 100,
mean( orientation_consistent_scm2 & significant_consistent_scm2 ) * 100
),digits=3)

orientation_consistent_tm = (((gamsub_tm[2,]<0) & (iso_sub$Estimate<0)) | ((gamsub_tm[2,]>0) & (iso_sub$Estimate>0)))
significant_consistent_tm = pvsub_tm[2,]<0.05
round(c(
sqrt(mean((iso_sub$Estimate - gamsub_tm[2,])^2)),
mean(significant_consistent_tm) * 100,
mean(orientation_consistent_tm) * 100,
mean( orientation_consistent_tm & significant_consistent_tm ) * 100
),digits=3)

orientation_consistent_cdmc = (((gamsub_cdmc[2,]<0) & (iso_sub$Estimate<0)) | ((gamsub_cdmc[2,]>0) & (iso_sub$Estimate>0)))
significant_consistent_cdmc = pvsub_cdmc[2,]<0.05
round(c(
sqrt(mean((iso_sub$Estimate - gamsub_cdmc[2,])^2)),
mean(significant_consistent_cdmc) * 100,
mean(orientation_consistent_cdmc) * 100,
mean( orientation_consistent_cdmc & significant_consistent_cdmc ) * 100
),digits=3)

orientation_consistent_tca = (((gamsub_tca[2,]<0) & (iso_sub$Estimate<0)) | ((gamsub_tca[2,]>0) & (iso_sub$Estimate>0)))
significant_consistent_tca = pvsub_tca[2,]<0.05
round(c(
sqrt(mean((iso_sub$Estimate - gamsub_tca[2,])^2)),
mean(significant_consistent_tca) * 100,
mean(orientation_consistent_tca) * 100,
mean( orientation_consistent_tca & significant_consistent_tca ) * 100
),digits=3)
