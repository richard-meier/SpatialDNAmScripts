
install.packages("gtools")
install.packages("R2OpenBUGS")

options(stringsAsFactors = FALSE)
library(gtools)
library(R2OpenBUGS)
library(MASS)

#############
# CONSTANTS #
################################################################################################################

num_subjects = 40 ### even sample sizes only !!!
delta.pow = 0.2
noise.z = 0.01
noise.x = 0.01

chrom = "chr1"
cluster_threshold1 = 3000
cluster_threshold2 = 30000

SIM_SEED = 11
CLUSTER_START = 1
MCMC_NUM = 15000
MCMC_BURN = 1000
MAX_CN = 10
MAX_SCN = 8
MAX_CPG_PER_SC = 44

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

set.seed(SIM_SEED)

source(paste0(inpath,"utility_scripts/utility_functions.R"))
source(paste0(inpath,"utility_scripts/two_arm_model_functions_sim.R"))
load(paste0(inpath,"data/DNA_methylation_data.RData"))

# create super cluster structures
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
rm(betas.mixtures)
rm(pheno.mixtures)
gc()


gammas = list()
TM_gamma_means = list()
TM_gamma_sds = list()
TM_psi_means = list()
TM_psi_sds = list()
TM_psi_medians = list()
TM_tau_means = list()
TM_tau_sds = list()
TM_tau_medians = list()
SCM2_gamma_means = list()
SCM2_gamma_sds = list()
SCM2_psi_means = list()
SCM2_psi_sds = list()
SCM2_psi_medians = list()
SCM2_tau_means = list()
SCM2_tau_sds = list()
SCM2_tau_medians = list()
SCM2_ogam_means = list()
SCM2_ogam_sds = list()
SCM2_oogam_means = list()
SCM2_oogam_sds = list()
SCM2_olambda_means = list()
SCM2_olambda_sds = list()
SCM2_oolambda_means = list()
SCM2_oolambda_sds = list()

results = data.frame(
	IT.ID=NA,MCMC.NUM=NA,MCMC.BURN=NA,MAX.SCN=NA,MAX.CN=NA,noise.z=NA,noise.x=NA,model=NA,fit.time=NA,num.clusters=NA,num.cpgs=NA,cluster.sizes=NA,
	super.cluster.span=NA,ncr.dev=NA,PD=NA,DIC=NA,RMSE=NA,PREDRMSE=NA
)
results = results[-1,]

for(rrr in CLUSTER_START:length(supercl)){ ##########################################################################################################################################################

	ctime = Sys.time()
	
	deltas = rep(1,6) * delta.pow
	nfacts = sample( c(-1,-1,0,0,1,1) )
	
	concentration_parameters = c(15.0727, 2.5392, 1.8439, 1.7934, 0.7404, 0.7240)
	cell_fractions = rdirichlet(n=num_subjects,alpha=concentration_parameters)
	cft = t(cell_fractions)
	cell_fractions_q = rdirichlet(n=num_subjects,alpha=concentration_parameters)
	cft_q = t(cell_fractions_q)
	colnames(cell_fractions) = colnames(cell_fractions_q) = c("Neutrophil", "Monocyte", "CD4T", "CD8T", "Bcell", "NK")
	num_ctypes = ncol(cell_fractions)
	
	clusters = supercl[[rrr]]
	max_clst_length = max(unlist(lapply(clusters,FUN=length)))
	
	pheno.sub = pheno.normal6[draw_random_mixture_candidates(pheno.normal6),]
	# pheno.sub order = c("Neutrophil", "NK", "Bcell", "CD4T", "CD8T","Monocyte")
	# order of concentration parameters = c("Neutrophil", "Monocyte", "CD4T", "CD8T", "Bcell", "NK")
	# manually reorder samples from most abundant to least abundant:
	pheno.sub = pheno.sub[c(1,6,4,5,3,2),]
	pheno.sub$CellType = as.character(pheno.sub$CellType)
	selected_samples = as.character(pheno.sub$Array)
	betas_sub = betas.normal6[,selected_samples ]
	
	
	GRPS = rep(0:1,each=num_subjects/2)
	X = array(NA,dim=c(num_subjects,length(clusters),max_clst_length))
	QPRED = array(NA,dim=c(num_subjects,length(clusters),max_clst_length))
	num_cpgs = c()
	
	cat(
		"processing super cluster",rrr,"/",length(supercl),"(",paste(sort(unlist(lapply(clusters,FUN=length))),collapse=","),") with",length(unlist(clusters)),
		"cpgs ; delta.pow =",delta.pow, "; sd(x,z) =",paste0("(",noise.x,",",noise.z,")"),
		"... \n"
	)
	
	allAN = anchr[anchr$pos %in% unlist(clusters),]
	allBE = betas_sub[row.names(betas_sub) %in% row.names(allAN),]
	ct_median = c()
	ct_mins = c()
	ct_maxs = c()
	cases = c()
	for(cti in 1:num_ctypes){
		ct_median[cti] = median(allBE[,cti])
		ct_mins[cti] = min(allBE[,cti])
		ct_maxs[cti] = max(allBE[,cti])
		
		if(ct_maxs[cti] < 0.6){
			cases[cti] = 1
		} else if(ct_mins[cti] > 0.4){
			cases[cti] = 2
		} else{
			if(ct_median[cti] < 0.5){
				cases[cti] = 1
			} else{ cases[cti] = 2 }
		}
	}
	
	tgamvec = c()
	for(clst in 1:length(clusters)){
		TNCPG = length(clusters[[clst]])
		for(cpgi in 1:TNCPG){
			for(cti in 1:num_ctypes){
				tgamvec = c(tgamvec,deltas[cti]*nfacts[cti])
			}
		}
	}
	true_gammas = matrix(tgamvec,nrow=6,byrow=FALSE)
	row.names(true_gammas) = colnames(cell_fractions)
	colnames(true_gammas) = anchr[anchr$pos %in% unlist(clusters),]$Name
	
	for(clst in 1:length(clusters)){
		TNCPG = length(clusters[[clst]])
		subAN = anchr[anchr$pos %in% clusters[[clst]],]
		subBE = betas_sub[row.names(betas_sub) %in% row.names(subAN),]
		subMixBE = matrix(NA,nrow=TNCPG,ncol=num_subjects)
		
		for(si in 1:num_subjects){
			tmpBE = matrix(subBE,nrow=TNCPG)
			for(cpgi in 1:TNCPG){ # CPG INDEX
				for(cti in 1:num_ctypes){ # CELL TYPE INDEX
					delta.gam = c(0,deltas[cti])
					if(cases[cti] == 2){
						delta.gam = c(-deltas[cti],0)
					}
					tmpBE[cpgi,cti] = tmpBE[cpgi,cti] + delta.gam[GRPS[si]+1]*nfacts[cti] + rnorm(n=1,mean=0,sd=noise.z)
				}
			}
			subMixBE[,si] = tmpBE %*% cft[,si] + rnorm(n=TNCPG,mean=0,sd=noise.x)
		}
		ctmp = t(subMixBE)
		X[1:nrow(ctmp),clst,1:ncol(ctmp)] = ctmp
		
		for(si in 1:num_subjects){
			tmpBE = matrix(subBE,nrow=TNCPG)
			for(cpgi in 1:TNCPG){ # CPG INDEX
				for(cti in 1:num_ctypes){ # CELL TYPE INDEX
					delta.gam = c(0,deltas[cti])
					if(cases[cti] == 2){
						delta.gam = c(-deltas[cti],0)
					}
					tmpBE[cpgi,cti] = tmpBE[cpgi,cti] + delta.gam[GRPS[si]+1]*nfacts[cti] + rnorm(n=1,mean=0,sd=noise.z)
				}
			}
			subMixBE[,si] = tmpBE %*% cft_q[,si] + rnorm(n=TNCPG,mean=0,sd=noise.x)
		}
		ctmp = t(subMixBE)
		QPRED[1:nrow(ctmp),clst,1:ncol(ctmp)] = ctmp
		
		num_cpgs[clst] = ncol(ctmp)
	}

	dats = list(
		X = X,
		Q = QPRED,
		GRP = GRPS,
		N = num_subjects,
		MCPG = num_cpgs,
		DEGCL = (num_cpgs <= 1)*1,
		KCF = num_ctypes,
		CFR = cell_fractions,
		QCFR = cell_fractions_q,
		MCL = length(clusters)
	)
	
	fit_times = list()
	inivals = generate.initial.values.sim(clusters=clusters, chains=1, KCF=6, N=num_subjects)
	
	
	ftime = Sys.time()
	fit_SCM2_cprior2 <- bugs(
		dats, inits = inivals, 
		parameters.to.save = c(
			"wrapper.mu","wrapper.gam","tau","psi","omu","oomu","ogam","oogam",
			"otau","ootau","olambda","oolambda", 
			"RMSE","PREDRMSE"
		), 
		model.file = SCM2_CP2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
		digits=8
	)
	fit_times[["SCM2.cprior2"]] = as.numeric(difftime(Sys.time(),ftime,units = "secs"))
	
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
	median_psi_SCM2_cp2 = matrix(
		fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="psi",5],
		nrow=6, byrow=TRUE
	)
	row.names(mean_psi_SCM2_cp2) = row.names(sd_psi_SCM2_cp2) = row.names(median_psi_SCM2_cp2) = colnames(cell_fractions)
	colnames(mean_psi_SCM2_cp2) = colnames(sd_psi_SCM2_cp2) = colnames(median_psi_SCM2_cp2) = anchr[anchr$pos %in% unlist(clusters),]$Name
	
	mean_olambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,7)=="olambda",1]
	sd_olambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,7)=="olambda",2]
	mean_oolambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,8)=="oolambda",1]
	sd_oolambda_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,8)=="oolambda",2]
	
	mean_tau_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="tau",1]
	sd_tau_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="tau",2]
	median_tau_SCM2_cp2 = fit_SCM2_cprior2$summary[substr(row.names(fit_SCM2_cprior2$summary),1,3)=="tau",5]


	ftime = Sys.time()
	fit_TM_cprior2 <- bugs(
		dats, inits = inivals, 
		parameters.to.save = c(
			"mu","gam","tau","psi","RMSE","PREDRMSE"
		), 
		model.file = TM_CP2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
		digits=8
	)
	fit_times[["TM.cprior2"]] = as.numeric(difftime(Sys.time(),ftime,units = "secs"))
	
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
	median_psi_TM_cp2 = matrix(
		fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="psi",5],
		nrow=6, byrow=TRUE
	)
	row.names(mean_psi_TM_cp2) = row.names(sd_psi_TM_cp2) = row.names(median_psi_TM_cp2) = colnames(cell_fractions)
	colnames(mean_psi_TM_cp2) = colnames(sd_psi_TM_cp2) = colnames(median_psi_TM_cp2) = anchr[anchr$pos %in% unlist(clusters),]$Name
	
	mean_tau_TM_cp2 = fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="tau",1]
	sd_tau_TM_cp2 = fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="tau",2]
	median_tau_TM_cp2 = fit_TM_cprior2$summary[substr(row.names(fit_TM_cprior2$summary),1,3)=="tau",5]
	
	
	i.res.tab = data.frame(
		IT.ID = rrr, MCMC.NUM=MCMC_NUM, MCMC.BURN=MCMC_BURN, MAX.SCN=MAX_SCN, MAX.CN=MAX_CN,
		noise.z=noise.z, noise.x=noise.x,
		model = c(
			"TM.cprior2", "SCM2.cprior2"
		),
		fit.time = c(
			fit_times[["TM.cprior2"]], fit_times[["SCM2.cprior2"]]
		),
		num.clusters = length(clusters),
		num.cpgs = length(unlist(clusters)),
		cluster.sizes = paste(sort(unlist(lapply(clusters,FUN=length))),collapse=","),
		super.cluster.span = calculate_super_cluster_span(clusters),
		ncr.dev = c(
			cross.rate(fit_TM_cprior2$sims.list[["deviance"]]),
			cross.rate(fit_SCM2_cprior2$sims.list[["deviance"]])
		),
		PD = c(
			fit_TM_cprior2$pD, fit_SCM2_cprior2$pD
		),
		DIC = c(
			fit_TM_cprior2$DIC, fit_SCM2_cprior2$DIC
		),
		RMSE = c(
			fit_TM_cprior2$summary[row.names(fit_TM_cprior2$summary) == "RMSE",1], 
			fit_SCM2_cprior2$summary[row.names(fit_SCM2_cprior2$summary) == "RMSE",1]
		),
		PREDRMSE = c(
			fit_TM_cprior2$summary[row.names(fit_TM_cprior2$summary) == "PREDRMSE",1], 
			fit_SCM2_cprior2$summary[row.names(fit_SCM2_cprior2$summary) == "PREDRMSE",1]
		)
	)
	results = rbind(results,i.res.tab)
	
	gammas[[rrr]] = true_gammas
	TM_gamma_means[[rrr]] = mean_gam_TM_cp2
	TM_gamma_sds[[rrr]] = sd_gam_TM_cp2
	TM_psi_means[[rrr]] = mean_psi_TM_cp2
	TM_psi_sds[[rrr]] = sd_psi_TM_cp2
	TM_psi_medians[[rrr]] = median_psi_TM_cp2
	TM_tau_means[[rrr]] = mean_tau_TM_cp2
	TM_tau_sds[[rrr]] = sd_tau_TM_cp2
	TM_tau_medians[[rrr]] = median_tau_TM_cp2
	SCM2_gamma_means[[rrr]] = mean_gam_SCM2_cp2
	SCM2_gamma_sds[[rrr]] = sd_gam_SCM2_cp2
	SCM2_psi_means[[rrr]] = mean_psi_SCM2_cp2
	SCM2_psi_sds[[rrr]] = sd_psi_SCM2_cp2
	SCM2_psi_medians[[rrr]] = median_psi_SCM2_cp2
	SCM2_tau_means[[rrr]] = mean_tau_SCM2_cp2
	SCM2_tau_sds[[rrr]] = sd_tau_SCM2_cp2
	SCM2_tau_medians[[rrr]] = median_tau_SCM2_cp2
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
		"simB2_cth1.",cluster_threshold1,"_cth2.",cluster_threshold2,"_",chrom,
		"_nsbj.",num_subjects,"_mcn.",MAX_CN,"_mscn.",MAX_SCN,"_delta.",delta.pow,
		"_noise.z.",noise.z,"_noise.x.",noise.x,"_seed.",SIM_SEED,"_clidx.",CLUSTER_START,
		"_mat.RDATA"
	)
	save(
		gammas, TM_gamma_means, TM_gamma_sds, TM_psi_means, TM_psi_sds, TM_psi_medians, TM_tau_means, TM_tau_sds, TM_tau_medians,
		SCM2_gamma_means, SCM2_gamma_sds, SCM2_psi_means, SCM2_psi_sds, SCM2_psi_medians, SCM2_tau_means, SCM2_tau_sds, SCM2_tau_medians,
		SCM2_ogam_means, SCM2_ogam_sds, SCM2_oogam_means, SCM2_oogam_sds,
		SCM2_olambda_means,SCM2_olambda_sds,SCM2_oolambda_means,SCM2_oolambda_sds,
		results,
		file=out_string
	)
	
	cat("   ... finished iteration after a "); print(Sys.time()-ctime)
	
	rm(
		dats, X, QPRED, cell_fractions, cell_fractions_q, cft, cft_q,
		fit_TM_cprior2,fit_SCM2_cprior2, pheno.sub, betas_sub, 
		subAN, subBE, subMixBE, tmpBE
	)
	gc()

} ##########################################################################################################################################################
