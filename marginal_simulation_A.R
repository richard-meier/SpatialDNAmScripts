options(stringsAsFactors = FALSE)

install.packages("gtools")
install.packages("R2OpenBUGS")

library(gtools)
library(R2OpenBUGS)
library(MASS)

#############
# CONSTANTS #
################################################################################################################

set.seed(11)
num_subjects = 20

noise.z = 0.05
noise.x = 0.05

chrom = "chr1"
cluster_threshold1 = 3000
cluster_threshold2 = 30000

MCMC_NUM = 15000
MCMC_BURN = 1000
MAX_CN = 10
MAX_SCN = 8
MAX_CPG_PER_SC = 44

inpath = ""
outpath = ""


#######################
# LOAD + EXECUTE ARGS #
################################################################################################################

args = commandArgs(TRUE)
for(i in 1:length(args)){
	eval(parse(text=args[[i]]))
}


################
# PROGRAM CODE #
################################################################################################################

source(paste0(inpath,"utility_scripts/utility_functions.R"))
source(paste0(inpath,"utility_scripts/marginal_model_specifications.R"))
source(paste0(inpath,"utility_scripts/initialize_models.R"))
load(paste0(inpath,"data/DNA_methylation_data.RData"))

anchr = annotEPIC
anchr = anchr[row.names(anchr) %in% row.names(betas.normal6),]
anchr = anchr[anchr$chr==chrom,]
anchr = anchr[order(anchr$pos),]

betas.normal6 = betas.normal6[row.names(betas.normal6) %in% row.names(anchr), ]

super_clusters = construct.cluster.structures(
	chr_positions=anchr$pos, cl_threshold1=cluster_threshold1, cl_threshold2=cluster_threshold2,
	max_cn=MAX_CN, max_scn=MAX_SCN, max_split_prop=0.7
)

supercl = super_clusters$sc.v3
cluster_group_length = unlist(lapply(supercl,FUN=length))
supercl = supercl[cluster_group_length > 1]

cpg_counts = unlist(lapply(supercl, FUN=function(x){return(length(unlist(x)))} ))
supercl = supercl[cpg_counts <= MAX_CPG_PER_SC] ### HARD SUBSET (save time and memory)



# release unused data from memory
rm(super_clusters)
rm(annotEPIC)
gc()



results = data.frame(
	IT.ID=NA, MCMC.NUM=NA, MCMC.BURN=NA, MAX.SCN=NA, MAX.CN=NA, noise.z=NA, noise.x=NA, 
	model=NA, num.clusters=NA, num.cpgs=NA, cluster.sizes=NA, super.cluster.span=NA, 
	ncr.dev=NA, PD=NA, DIC=NA, RMSE=NA, PREDRMSE=NA
)
results = results[-1,]

for(rrr in 1:length(supercl)){ ######################################################################################################

ctime = Sys.time()

cell_fractions = rdirichlet(n=num_subjects,alpha=rep(2,6))
cft = t(cell_fractions)
num_ctypes = ncol(cell_fractions)

cell_fractions_q = rdirichlet(n=num_subjects,alpha=rep(2,6))
cft_q = t(cell_fractions_q)

clusters = supercl[[rrr]]
max_clst_length = max(unlist(lapply(clusters,FUN=length)))

pheno.sub = pheno.normal6[draw_random_mixture_candidates(pheno.normal6),]
selected_samples = as.character(pheno.sub$Array)
betas_sub = betas.normal6[,selected_samples ]

X = array(NA,dim=c(num_subjects,length(clusters),max_clst_length))
QPRED = array(NA,dim=c(num_subjects,length(clusters),max_clst_length))
num_cpgs = c()



cat("processing super cluster",rrr,"/",length(supercl),"(",paste(sort(unlist(lapply(clusters,FUN=length))),collapse=","),") with",length(unlist(clusters)),"cpgs ...\n")
	
for(clst in 1:length(clusters)){
	TNCPG = length(clusters[[clst]])
	subAN = anchr[anchr$pos %in% clusters[[clst]],]
	subBE = betas_sub[row.names(betas_sub) %in% row.names(subAN),]
	subMixBE = matrix(NA,nrow=TNCPG,ncol=num_subjects)
	
	for(si in 1:num_subjects){
		tmpBE = matrix(subBE,nrow=TNCPG)
		for(sss1 in 1:TNCPG){ # CPG INDEX
			for(sss2 in 1:num_ctypes){ # CELL TYPE INDEX
				tmpBE[sss1,sss2] = tmpBE[sss1,sss2] + rnorm(n=1,mean=0,sd=noise.z)
			}
		}
		subMixBE[,si] = tmpBE %*% cft[,si] + rnorm(n=1,mean=0,sd=noise.x)
	}
	ctmp = t(subMixBE)
	X[1:nrow(ctmp),clst,1:ncol(ctmp)] = ctmp
	
	for(si in 1:num_subjects){
		tmpBE = matrix(subBE,nrow=TNCPG)
		for(sss1 in 1:TNCPG){ # CPG INDEX
			for(sss2 in 1:num_ctypes){ # CELL TYPE INDEX
				tmpBE[sss1,sss2] = tmpBE[sss1,sss2] + rnorm(n=1,mean=0,sd=noise.z)
			}
		}
		subMixBE[,si] = tmpBE %*% cft_q[,si] + rnorm(n=1,mean=0,sd=noise.x)
	}
	ctmp = t(subMixBE)
	QPRED[1:nrow(ctmp),clst,1:ncol(ctmp)] = ctmp
	
	num_cpgs[clst] = ncol(ctmp)
}

dats = list(
	X = X,
	Q = QPRED,
	N = num_subjects,
	MCPG = num_cpgs,
	DEGCL = (num_cpgs <= 1)*1,
	KCF = num_ctypes,
	CFR = cell_fractions,
	QCFR = cell_fractions_q,
	MCL = length(clusters)
)

inivals = generate.initial.values.berry.m1(clusters=clusters, chains = 1, KCF=6, N=num_subjects)

fit_berry_m1 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","tau","psi","omu","oomu","ooomu","otau","ootau","oootau",
		"RMSE","PREDRMSE"
	), 
	model.file = marginal_berry_like_model_m1, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_berry_m1_cprior <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","tau","psi","omu","oomu","ooomu","otau","ootau","oootau",
		"RMSE","PREDRMSE"
	), 
	model.file = marginal_berry_like_model_m1_custom_prior, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_berry_m1_cprior2 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","tau","psi","omu","oomu","ooomu","otau","ootau","oootau",
		"RMSE","PREDRMSE"
	), 
	model.file = marginal_berry_like_model_m1_custom_prior2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_berry_m2 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","tau","psi","omu","oomu","otau","ootau",
		"RMSE","PREDRMSE"
	), 
	model.file = marginal_berry_like_model_m2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_berry_m2_cprior <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","tau","psi","omu","oomu","otau","ootau",
		"RMSE","PREDRMSE"
	), 
	model.file = marginal_berry_like_model_m2_custom_prior, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_berry_m2_cprior2 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","tau","psi","omu","oomu","otau","ootau",
		"RMSE","PREDRMSE"
	), 
	model.file = marginal_berry_like_model_m2_custom_prior2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_berry_m3 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","tau","psi","oomu","ootau",
		"RMSE","PREDRMSE"
	), 
	model.file = marginal_berry_like_model_m3, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_berry_m3_cprior <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","tau","psi","oomu","ootau",
		"RMSE","PREDRMSE"
	), 
	model.file = marginal_berry_like_model_m3_custom_prior, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_berry_m3_cprior2 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","tau","psi","oomu","ootau",
		"RMSE","PREDRMSE"
	), 
	model.file = marginal_berry_like_model_m3_custom_prior2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_tca <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"mu","tau","psi","RMSE","PREDRMSE"
	), 
	model.file = marginal_tca_like_model_flip, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_tca_cprior <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"mu","tau","psi","RMSE","PREDRMSE"
	), 
	model.file = marginal_tca_like_model_flip_custom_prior, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

fit_tca_cprior2 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"mu","tau","psi","RMSE","PREDRMSE"
	), 
	model.file = marginal_tca_like_model_flip_custom_prior2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)

i.res.tab = data.frame(
	IT.ID = rrr, MCMC.NUM=MCMC_NUM, MCMC.BURN=MCMC_BURN, MAX.SCN=MAX_SCN, MAX.CN=MAX_CN,
	noise.z=noise.z, noise.x=noise.x, 
	model = c(
		"TCA", "TCA.cprior1", "TCA.cprior2",
		"Berry.M1", "Berry.M1.cprior1", "Berry.M1.cprior2",
		"Berry.M2", "Berry.M2.cprior1", "Berry.M2.cprior2",
		"Berry.M3", "Berry.M3.cprior1", "Berry.M3.cprior2"
	),
	num.clusters = length(clusters),
	num.cpgs = length(unlist(clusters)),
	cluster.sizes = paste(sort(unlist(lapply(clusters,FUN=length))),collapse=","),
	super.cluster.span = calculate_super_cluster_span(clusters),
	ncr.dev = c(
		cross.rate(fit_tca$sims.list[["deviance"]]),
		cross.rate(fit_tca_cprior$sims.list[["deviance"]]),
		cross.rate(fit_tca_cprior2$sims.list[["deviance"]]),
		cross.rate(fit_berry_m1$sims.list[["deviance"]]),
		cross.rate(fit_berry_m1_cprior$sims.list[["deviance"]]),
		cross.rate(fit_berry_m1_cprior2$sims.list[["deviance"]]),
		cross.rate(fit_berry_m2$sims.list[["deviance"]]),
		cross.rate(fit_berry_m2_cprior$sims.list[["deviance"]]),
		cross.rate(fit_berry_m2_cprior2$sims.list[["deviance"]]),
		cross.rate(fit_berry_m3$sims.list[["deviance"]]),
		cross.rate(fit_berry_m3_cprior$sims.list[["deviance"]]),
		cross.rate(fit_berry_m3_cprior2$sims.list[["deviance"]])
	),
	PD = c(
		fit_tca$pD, fit_tca_cprior$pD, fit_tca_cprior2$pD, 
		fit_berry_m1$pD, fit_berry_m1_cprior$pD, fit_berry_m1_cprior2$pD, 
		fit_berry_m2$pD, fit_berry_m2_cprior$pD, fit_berry_m2_cprior2$pD, 
		fit_berry_m3$pD, fit_berry_m3_cprior$pD, fit_berry_m3_cprior2$pD
	),
	DIC = c(
		fit_tca$DIC, fit_tca_cprior$DIC, fit_tca_cprior2$DIC, 
		fit_berry_m1$DIC, fit_berry_m1_cprior$DIC, fit_berry_m1_cprior2$DIC, 
		fit_berry_m2$DIC, fit_berry_m2_cprior$DIC, fit_berry_m2_cprior2$DIC, 
		fit_berry_m3$DIC, fit_berry_m3_cprior$DIC, fit_berry_m3_cprior2$DIC
	),
	RMSE = c(
		fit_tca$summary[row.names(fit_tca$summary) == "RMSE",1], 
		fit_tca_cprior$summary[row.names(fit_tca_cprior$summary) == "RMSE",1], 
		fit_tca_cprior2$summary[row.names(fit_tca_cprior2$summary) == "RMSE",1], 
		fit_berry_m1$summary[row.names(fit_berry_m1$summary) == "RMSE",1],
		fit_berry_m1_cprior$summary[row.names(fit_berry_m1_cprior$summary) == "RMSE",1],
		fit_berry_m1_cprior2$summary[row.names(fit_berry_m1_cprior2$summary) == "RMSE",1],
		fit_berry_m2$summary[row.names(fit_berry_m2$summary) == "RMSE",1],
		fit_berry_m2_cprior$summary[row.names(fit_berry_m2_cprior$summary) == "RMSE",1],
		fit_berry_m2_cprior2$summary[row.names(fit_berry_m2_cprior2$summary) == "RMSE",1],
		fit_berry_m3$summary[row.names(fit_berry_m3$summary) == "RMSE",1],
		fit_berry_m3_cprior$summary[row.names(fit_berry_m3_cprior$summary) == "RMSE",1],
		fit_berry_m3_cprior2$summary[row.names(fit_berry_m3_cprior2$summary) == "RMSE",1]
	),
	PREDRMSE = c(
		fit_tca$summary[row.names(fit_tca$summary) == "PREDRMSE",1], 
		fit_tca_cprior$summary[row.names(fit_tca_cprior$summary) == "PREDRMSE",1], 
		fit_tca_cprior2$summary[row.names(fit_tca_cprior2$summary) == "PREDRMSE",1], 
		fit_berry_m1$summary[row.names(fit_berry_m1$summary) == "PREDRMSE",1],
		fit_berry_m1_cprior$summary[row.names(fit_berry_m1_cprior$summary) == "PREDRMSE",1],
		fit_berry_m1_cprior2$summary[row.names(fit_berry_m1_cprior2$summary) == "PREDRMSE",1],
		fit_berry_m2$summary[row.names(fit_berry_m2$summary) == "PREDRMSE",1],
		fit_berry_m2_cprior$summary[row.names(fit_berry_m2_cprior$summary) == "PREDRMSE",1],
		fit_berry_m2_cprior2$summary[row.names(fit_berry_m2_cprior2$summary) == "PREDRMSE",1],
		fit_berry_m3$summary[row.names(fit_berry_m3$summary) == "PREDRMSE",1],
		fit_berry_m3_cprior$summary[row.names(fit_berry_m3_cprior$summary) == "PREDRMSE",1],
		fit_berry_m3_cprior2$summary[row.names(fit_berry_m3_cprior2$summary) == "PREDRMSE",1]
	)
)

results = rbind(results,i.res.tab)

out_string = paste0(
	outpath,
	"marginal_simA_cth1.",cluster_threshold1,"_cth2.",cluster_threshold2,"_",chrom,
	"_nsbj.",num_subjects,"_mcn.",MAX_CN,"_mscn.",MAX_SCN,"_noise.z",noise.z,"_noise.x",noise.x,"_table.csv"
)
write.csv(results,file=out_string,row.names=FALSE)

cat("   ... finished iteration after a "); print(Sys.time()-ctime)

rm(
	dats, X, QPRED, cell_fractions, cell_fractions_q, cft, cft_q,
	fit_tca, fit_tca_cprior, fit_tca_cprior2, 
	fit_berry_m1, fit_berry_m1_cprior, fit_berry_m1_cprior2, 
	fit_berry_m2, fit_berry_m2_cprior, fit_berry_m2_cprior2, 
	fit_berry_m3, fit_berry_m3_cprior, fit_berry_m3_cprior2,
	betas_sub, subAN, subBE, subMixBE, tmpBE
)
gc()


} ##########################################################################################################

