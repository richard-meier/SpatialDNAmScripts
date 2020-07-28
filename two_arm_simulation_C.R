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
num_subjects = 40 ### even sample sizes only !!!
delta.pow = 0.3
noise.z = 0.1
noise.x = 0.1
noise.c = 0.001
cor.c.w = 0.5 		# within cluster
cor.c.b = cor.c.w/2 	# between cluster

chrom = "chr1"
cluster_threshold1 = 3000
cluster_threshold2 = 30000

tquants = seq(from=0.2,to=0.001,by=-0.001) #(20:1) * 0.01

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

var.cluster = (noise.c)^2


################
# PROGRAM CODE #
################################################################################################################

source(paste0(inpath,"utility_scripts/utility_functions.R"))
source(paste0(inpath,"utility_scripts/model_specifications.R"))
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
	IT.ID=NA,MCMC.NUM=NA,MCMC.BURN=NA,MAX.SCN=NA,MAX.CN=NA,noise.z=NA,noise.x=NA,noise.c=NA,cor.c.w=NA,cor.c.b=NA,model=NA,fit.time=NA,num.clusters=NA,num.cpgs=NA,cluster.sizes=NA,
super.cluster.span=NA,ncr.dev=NA,PD=NA,DIC=NA,RMSE=NA,PREDRMSE=NA,pw.q.200=NA,pw.q.199=NA,pw.q.198=NA,pw.q.197=NA,pw.q.196=NA,pw.q.195=NA,pw.q.194=NA,pw.q.193=NA,pw.q.192=NA,
pw.q.191=NA,pw.q.190=NA,pw.q.189=NA,pw.q.188=NA,pw.q.187=NA,pw.q.186=NA,pw.q.185=NA,pw.q.184=NA,pw.q.183=NA,pw.q.182=NA,pw.q.181=NA,pw.q.180=NA,pw.q.179=NA,pw.q.178=NA,pw.q.177=NA,
pw.q.176=NA,pw.q.175=NA,pw.q.174=NA,pw.q.173=NA,pw.q.172=NA,pw.q.171=NA,pw.q.170=NA,pw.q.169=NA,pw.q.168=NA,pw.q.167=NA,pw.q.166=NA,pw.q.165=NA,pw.q.164=NA,pw.q.163=NA,pw.q.162=NA,
pw.q.161=NA,pw.q.160=NA,pw.q.159=NA,pw.q.158=NA,pw.q.157=NA,pw.q.156=NA,pw.q.155=NA,pw.q.154=NA,pw.q.153=NA,pw.q.152=NA,pw.q.151=NA,pw.q.150=NA,pw.q.149=NA,pw.q.148=NA,pw.q.147=NA,
pw.q.146=NA,pw.q.145=NA,pw.q.144=NA,pw.q.143=NA,pw.q.142=NA,pw.q.141=NA,pw.q.140=NA,pw.q.139=NA,pw.q.138=NA,pw.q.137=NA,pw.q.136=NA,pw.q.135=NA,pw.q.134=NA,pw.q.133=NA,pw.q.132=NA,
pw.q.131=NA,pw.q.130=NA,pw.q.129=NA,pw.q.128=NA,pw.q.127=NA,pw.q.126=NA,pw.q.125=NA,pw.q.124=NA,pw.q.123=NA,pw.q.122=NA,pw.q.121=NA,pw.q.120=NA,pw.q.119=NA,pw.q.118=NA,pw.q.117=NA,
pw.q.116=NA,pw.q.115=NA,pw.q.114=NA,pw.q.113=NA,pw.q.112=NA,pw.q.111=NA,pw.q.110=NA,pw.q.109=NA,pw.q.108=NA,pw.q.107=NA,pw.q.106=NA,pw.q.105=NA,pw.q.104=NA,pw.q.103=NA,pw.q.102=NA,
pw.q.101=NA,pw.q.100=NA,pw.q.099=NA,pw.q.098=NA,pw.q.097=NA,pw.q.096=NA,pw.q.095=NA,pw.q.094=NA,pw.q.093=NA,pw.q.092=NA,pw.q.091=NA,pw.q.090=NA,pw.q.089=NA,pw.q.088=NA,pw.q.087=NA,
pw.q.086=NA,pw.q.085=NA,pw.q.084=NA,pw.q.083=NA,pw.q.082=NA,pw.q.081=NA,pw.q.080=NA,pw.q.079=NA,pw.q.078=NA,pw.q.077=NA,pw.q.076=NA,pw.q.075=NA,pw.q.074=NA,pw.q.073=NA,pw.q.072=NA,
pw.q.071=NA,pw.q.070=NA,pw.q.069=NA,pw.q.068=NA,pw.q.067=NA,pw.q.066=NA,pw.q.065=NA,pw.q.064=NA,pw.q.063=NA,pw.q.062=NA,pw.q.061=NA,pw.q.060=NA,pw.q.059=NA,pw.q.058=NA,pw.q.057=NA,
pw.q.056=NA,pw.q.055=NA,pw.q.054=NA,pw.q.053=NA,pw.q.052=NA,pw.q.051=NA,pw.q.050=NA,pw.q.049=NA,pw.q.048=NA,pw.q.047=NA,pw.q.046=NA,pw.q.045=NA,pw.q.044=NA,pw.q.043=NA,pw.q.042=NA,
pw.q.041=NA,pw.q.040=NA,pw.q.039=NA,pw.q.038=NA,pw.q.037=NA,pw.q.036=NA,pw.q.035=NA,pw.q.034=NA,pw.q.033=NA,pw.q.032=NA,pw.q.031=NA,pw.q.030=NA,pw.q.029=NA,pw.q.028=NA,pw.q.027=NA,
pw.q.026=NA,pw.q.025=NA,pw.q.024=NA,pw.q.023=NA,pw.q.022=NA,pw.q.021=NA,pw.q.020=NA,pw.q.019=NA,pw.q.018=NA,pw.q.017=NA,pw.q.016=NA,pw.q.015=NA,pw.q.014=NA,pw.q.013=NA,pw.q.012=NA,
pw.q.011=NA,pw.q.010=NA,pw.q.009=NA,pw.q.008=NA,pw.q.007=NA,pw.q.006=NA,pw.q.005=NA,pw.q.004=NA,pw.q.003=NA,pw.q.002=NA,pw.q.001=NA,cnt.q.200=NA,cnt.q.199=NA,cnt.q.198=NA,cnt.q.197=NA,
cnt.q.196=NA,cnt.q.195=NA,cnt.q.194=NA,cnt.q.193=NA,cnt.q.192=NA,cnt.q.191=NA,cnt.q.190=NA,cnt.q.189=NA,cnt.q.188=NA,cnt.q.187=NA,cnt.q.186=NA,cnt.q.185=NA,cnt.q.184=NA,cnt.q.183=NA,cnt.q.182=NA,
cnt.q.181=NA,cnt.q.180=NA,cnt.q.179=NA,cnt.q.178=NA,cnt.q.177=NA,cnt.q.176=NA,cnt.q.175=NA,cnt.q.174=NA,cnt.q.173=NA,cnt.q.172=NA,cnt.q.171=NA,cnt.q.170=NA,cnt.q.169=NA,cnt.q.168=NA,cnt.q.167=NA,
cnt.q.166=NA,cnt.q.165=NA,cnt.q.164=NA,cnt.q.163=NA,cnt.q.162=NA,cnt.q.161=NA,cnt.q.160=NA,cnt.q.159=NA,cnt.q.158=NA,cnt.q.157=NA,cnt.q.156=NA,cnt.q.155=NA,cnt.q.154=NA,cnt.q.153=NA,cnt.q.152=NA,
cnt.q.151=NA,cnt.q.150=NA,cnt.q.149=NA,cnt.q.148=NA,cnt.q.147=NA,cnt.q.146=NA,cnt.q.145=NA,cnt.q.144=NA,cnt.q.143=NA,cnt.q.142=NA,cnt.q.141=NA,cnt.q.140=NA,cnt.q.139=NA,cnt.q.138=NA,cnt.q.137=NA,
cnt.q.136=NA,cnt.q.135=NA,cnt.q.134=NA,cnt.q.133=NA,cnt.q.132=NA,cnt.q.131=NA,cnt.q.130=NA,cnt.q.129=NA,cnt.q.128=NA,cnt.q.127=NA,cnt.q.126=NA,cnt.q.125=NA,cnt.q.124=NA,cnt.q.123=NA,cnt.q.122=NA,
cnt.q.121=NA,cnt.q.120=NA,cnt.q.119=NA,cnt.q.118=NA,cnt.q.117=NA,cnt.q.116=NA,cnt.q.115=NA,cnt.q.114=NA,cnt.q.113=NA,cnt.q.112=NA,cnt.q.111=NA,cnt.q.110=NA,cnt.q.109=NA,cnt.q.108=NA,cnt.q.107=NA,
cnt.q.106=NA,cnt.q.105=NA,cnt.q.104=NA,cnt.q.103=NA,cnt.q.102=NA,cnt.q.101=NA,cnt.q.100=NA,cnt.q.099=NA,cnt.q.098=NA,cnt.q.097=NA,cnt.q.096=NA,cnt.q.095=NA,cnt.q.094=NA,cnt.q.093=NA,cnt.q.092=NA,
cnt.q.091=NA,cnt.q.090=NA,cnt.q.089=NA,cnt.q.088=NA,cnt.q.087=NA,cnt.q.086=NA,cnt.q.085=NA,cnt.q.084=NA,cnt.q.083=NA,cnt.q.082=NA,cnt.q.081=NA,cnt.q.080=NA,cnt.q.079=NA,cnt.q.078=NA,cnt.q.077=NA,
cnt.q.076=NA,cnt.q.075=NA,cnt.q.074=NA,cnt.q.073=NA,cnt.q.072=NA,cnt.q.071=NA,cnt.q.070=NA,cnt.q.069=NA,cnt.q.068=NA,cnt.q.067=NA,cnt.q.066=NA,cnt.q.065=NA,cnt.q.064=NA,cnt.q.063=NA,cnt.q.062=NA,
cnt.q.061=NA,cnt.q.060=NA,cnt.q.059=NA,cnt.q.058=NA,cnt.q.057=NA,cnt.q.056=NA,cnt.q.055=NA,cnt.q.054=NA,cnt.q.053=NA,cnt.q.052=NA,cnt.q.051=NA,cnt.q.050=NA,cnt.q.049=NA,cnt.q.048=NA,cnt.q.047=NA,
cnt.q.046=NA,cnt.q.045=NA,cnt.q.044=NA,cnt.q.043=NA,cnt.q.042=NA,cnt.q.041=NA,cnt.q.040=NA,cnt.q.039=NA,cnt.q.038=NA,cnt.q.037=NA,cnt.q.036=NA,cnt.q.035=NA,cnt.q.034=NA,cnt.q.033=NA,cnt.q.032=NA,
cnt.q.031=NA,cnt.q.030=NA,cnt.q.029=NA,cnt.q.028=NA,cnt.q.027=NA,cnt.q.026=NA,cnt.q.025=NA,cnt.q.024=NA,cnt.q.023=NA,cnt.q.022=NA,cnt.q.021=NA,cnt.q.020=NA,cnt.q.019=NA,cnt.q.018=NA,cnt.q.017=NA,
cnt.q.016=NA,cnt.q.015=NA,cnt.q.014=NA,cnt.q.013=NA,cnt.q.012=NA,cnt.q.011=NA,cnt.q.010=NA,cnt.q.009=NA,cnt.q.008=NA,cnt.q.007=NA,cnt.q.006=NA,cnt.q.005=NA,cnt.q.004=NA,cnt.q.003=NA,cnt.q.002=NA,
cnt.q.001=NA,CT1.pw.q.200=NA,CT2.pw.q.200=NA,CT3.pw.q.200=NA,CT4.pw.q.200=NA,CT5.pw.q.200=NA,CT6.pw.q.200=NA,CT1.pw.q.199=NA,CT2.pw.q.199=NA,CT3.pw.q.199=NA,CT4.pw.q.199=NA,CT5.pw.q.199=NA,CT6.pw.q.199=NA,CT1.pw.q.198=NA,CT2.pw.q.198=NA,
CT3.pw.q.198=NA,CT4.pw.q.198=NA,CT5.pw.q.198=NA,CT6.pw.q.198=NA,CT1.pw.q.197=NA,CT2.pw.q.197=NA,CT3.pw.q.197=NA,CT4.pw.q.197=NA,CT5.pw.q.197=NA,CT6.pw.q.197=NA,CT1.pw.q.196=NA,CT2.pw.q.196=NA,CT3.pw.q.196=NA,CT4.pw.q.196=NA,CT5.pw.q.196=NA,
CT6.pw.q.196=NA,CT1.pw.q.195=NA,CT2.pw.q.195=NA,CT3.pw.q.195=NA,CT4.pw.q.195=NA,CT5.pw.q.195=NA,CT6.pw.q.195=NA,CT1.pw.q.194=NA,CT2.pw.q.194=NA,CT3.pw.q.194=NA,CT4.pw.q.194=NA,CT5.pw.q.194=NA,CT6.pw.q.194=NA,CT1.pw.q.193=NA,CT2.pw.q.193=NA,
CT3.pw.q.193=NA,CT4.pw.q.193=NA,CT5.pw.q.193=NA,CT6.pw.q.193=NA,CT1.pw.q.192=NA,CT2.pw.q.192=NA,CT3.pw.q.192=NA,CT4.pw.q.192=NA,CT5.pw.q.192=NA,CT6.pw.q.192=NA,CT1.pw.q.191=NA,CT2.pw.q.191=NA,CT3.pw.q.191=NA,CT4.pw.q.191=NA,CT5.pw.q.191=NA,
CT6.pw.q.191=NA,CT1.pw.q.190=NA,CT2.pw.q.190=NA,CT3.pw.q.190=NA,CT4.pw.q.190=NA,CT5.pw.q.190=NA,CT6.pw.q.190=NA,CT1.pw.q.189=NA,CT2.pw.q.189=NA,CT3.pw.q.189=NA,CT4.pw.q.189=NA,CT5.pw.q.189=NA,CT6.pw.q.189=NA,CT1.pw.q.188=NA,CT2.pw.q.188=NA,
CT3.pw.q.188=NA,CT4.pw.q.188=NA,CT5.pw.q.188=NA,CT6.pw.q.188=NA,CT1.pw.q.187=NA,CT2.pw.q.187=NA,CT3.pw.q.187=NA,CT4.pw.q.187=NA,CT5.pw.q.187=NA,CT6.pw.q.187=NA,CT1.pw.q.186=NA,CT2.pw.q.186=NA,CT3.pw.q.186=NA,CT4.pw.q.186=NA,CT5.pw.q.186=NA,
CT6.pw.q.186=NA,CT1.pw.q.185=NA,CT2.pw.q.185=NA,CT3.pw.q.185=NA,CT4.pw.q.185=NA,CT5.pw.q.185=NA,CT6.pw.q.185=NA,CT1.pw.q.184=NA,CT2.pw.q.184=NA,CT3.pw.q.184=NA,CT4.pw.q.184=NA,CT5.pw.q.184=NA,CT6.pw.q.184=NA,CT1.pw.q.183=NA,CT2.pw.q.183=NA,
CT3.pw.q.183=NA,CT4.pw.q.183=NA,CT5.pw.q.183=NA,CT6.pw.q.183=NA,CT1.pw.q.182=NA,CT2.pw.q.182=NA,CT3.pw.q.182=NA,CT4.pw.q.182=NA,CT5.pw.q.182=NA,CT6.pw.q.182=NA,CT1.pw.q.181=NA,CT2.pw.q.181=NA,CT3.pw.q.181=NA,CT4.pw.q.181=NA,CT5.pw.q.181=NA,
CT6.pw.q.181=NA,CT1.pw.q.180=NA,CT2.pw.q.180=NA,CT3.pw.q.180=NA,CT4.pw.q.180=NA,CT5.pw.q.180=NA,CT6.pw.q.180=NA,CT1.pw.q.179=NA,CT2.pw.q.179=NA,CT3.pw.q.179=NA,CT4.pw.q.179=NA,CT5.pw.q.179=NA,CT6.pw.q.179=NA,CT1.pw.q.178=NA,CT2.pw.q.178=NA,
CT3.pw.q.178=NA,CT4.pw.q.178=NA,CT5.pw.q.178=NA,CT6.pw.q.178=NA,CT1.pw.q.177=NA,CT2.pw.q.177=NA,CT3.pw.q.177=NA,CT4.pw.q.177=NA,CT5.pw.q.177=NA,CT6.pw.q.177=NA,CT1.pw.q.176=NA,CT2.pw.q.176=NA,CT3.pw.q.176=NA,CT4.pw.q.176=NA,CT5.pw.q.176=NA,
CT6.pw.q.176=NA,CT1.pw.q.175=NA,CT2.pw.q.175=NA,CT3.pw.q.175=NA,CT4.pw.q.175=NA,CT5.pw.q.175=NA,CT6.pw.q.175=NA,CT1.pw.q.174=NA,CT2.pw.q.174=NA,CT3.pw.q.174=NA,CT4.pw.q.174=NA,CT5.pw.q.174=NA,CT6.pw.q.174=NA,CT1.pw.q.173=NA,CT2.pw.q.173=NA,
CT3.pw.q.173=NA,CT4.pw.q.173=NA,CT5.pw.q.173=NA,CT6.pw.q.173=NA,CT1.pw.q.172=NA,CT2.pw.q.172=NA,CT3.pw.q.172=NA,CT4.pw.q.172=NA,CT5.pw.q.172=NA,CT6.pw.q.172=NA,CT1.pw.q.171=NA,CT2.pw.q.171=NA,CT3.pw.q.171=NA,CT4.pw.q.171=NA,CT5.pw.q.171=NA,
CT6.pw.q.171=NA,CT1.pw.q.170=NA,CT2.pw.q.170=NA,CT3.pw.q.170=NA,CT4.pw.q.170=NA,CT5.pw.q.170=NA,CT6.pw.q.170=NA,CT1.pw.q.169=NA,CT2.pw.q.169=NA,CT3.pw.q.169=NA,CT4.pw.q.169=NA,CT5.pw.q.169=NA,CT6.pw.q.169=NA,CT1.pw.q.168=NA,CT2.pw.q.168=NA,
CT3.pw.q.168=NA,CT4.pw.q.168=NA,CT5.pw.q.168=NA,CT6.pw.q.168=NA,CT1.pw.q.167=NA,CT2.pw.q.167=NA,CT3.pw.q.167=NA,CT4.pw.q.167=NA,CT5.pw.q.167=NA,CT6.pw.q.167=NA,CT1.pw.q.166=NA,CT2.pw.q.166=NA,CT3.pw.q.166=NA,CT4.pw.q.166=NA,CT5.pw.q.166=NA,
CT6.pw.q.166=NA,CT1.pw.q.165=NA,CT2.pw.q.165=NA,CT3.pw.q.165=NA,CT4.pw.q.165=NA,CT5.pw.q.165=NA,CT6.pw.q.165=NA,CT1.pw.q.164=NA,CT2.pw.q.164=NA,CT3.pw.q.164=NA,CT4.pw.q.164=NA,CT5.pw.q.164=NA,CT6.pw.q.164=NA,CT1.pw.q.163=NA,CT2.pw.q.163=NA,
CT3.pw.q.163=NA,CT4.pw.q.163=NA,CT5.pw.q.163=NA,CT6.pw.q.163=NA,CT1.pw.q.162=NA,CT2.pw.q.162=NA,CT3.pw.q.162=NA,CT4.pw.q.162=NA,CT5.pw.q.162=NA,CT6.pw.q.162=NA,CT1.pw.q.161=NA,CT2.pw.q.161=NA,CT3.pw.q.161=NA,CT4.pw.q.161=NA,CT5.pw.q.161=NA,
CT6.pw.q.161=NA,CT1.pw.q.160=NA,CT2.pw.q.160=NA,CT3.pw.q.160=NA,CT4.pw.q.160=NA,CT5.pw.q.160=NA,CT6.pw.q.160=NA,CT1.pw.q.159=NA,CT2.pw.q.159=NA,CT3.pw.q.159=NA,CT4.pw.q.159=NA,CT5.pw.q.159=NA,CT6.pw.q.159=NA,CT1.pw.q.158=NA,CT2.pw.q.158=NA,
CT3.pw.q.158=NA,CT4.pw.q.158=NA,CT5.pw.q.158=NA,CT6.pw.q.158=NA,CT1.pw.q.157=NA,CT2.pw.q.157=NA,CT3.pw.q.157=NA,CT4.pw.q.157=NA,CT5.pw.q.157=NA,CT6.pw.q.157=NA,CT1.pw.q.156=NA,CT2.pw.q.156=NA,CT3.pw.q.156=NA,CT4.pw.q.156=NA,CT5.pw.q.156=NA,
CT6.pw.q.156=NA,CT1.pw.q.155=NA,CT2.pw.q.155=NA,CT3.pw.q.155=NA,CT4.pw.q.155=NA,CT5.pw.q.155=NA,CT6.pw.q.155=NA,CT1.pw.q.154=NA,CT2.pw.q.154=NA,CT3.pw.q.154=NA,CT4.pw.q.154=NA,CT5.pw.q.154=NA,CT6.pw.q.154=NA,CT1.pw.q.153=NA,CT2.pw.q.153=NA,
CT3.pw.q.153=NA,CT4.pw.q.153=NA,CT5.pw.q.153=NA,CT6.pw.q.153=NA,CT1.pw.q.152=NA,CT2.pw.q.152=NA,CT3.pw.q.152=NA,CT4.pw.q.152=NA,CT5.pw.q.152=NA,CT6.pw.q.152=NA,CT1.pw.q.151=NA,CT2.pw.q.151=NA,CT3.pw.q.151=NA,CT4.pw.q.151=NA,CT5.pw.q.151=NA,
CT6.pw.q.151=NA,CT1.pw.q.150=NA,CT2.pw.q.150=NA,CT3.pw.q.150=NA,CT4.pw.q.150=NA,CT5.pw.q.150=NA,CT6.pw.q.150=NA,CT1.pw.q.149=NA,CT2.pw.q.149=NA,CT3.pw.q.149=NA,CT4.pw.q.149=NA,CT5.pw.q.149=NA,CT6.pw.q.149=NA,CT1.pw.q.148=NA,CT2.pw.q.148=NA,
CT3.pw.q.148=NA,CT4.pw.q.148=NA,CT5.pw.q.148=NA,CT6.pw.q.148=NA,CT1.pw.q.147=NA,CT2.pw.q.147=NA,CT3.pw.q.147=NA,CT4.pw.q.147=NA,CT5.pw.q.147=NA,CT6.pw.q.147=NA,CT1.pw.q.146=NA,CT2.pw.q.146=NA,CT3.pw.q.146=NA,CT4.pw.q.146=NA,CT5.pw.q.146=NA,
CT6.pw.q.146=NA,CT1.pw.q.145=NA,CT2.pw.q.145=NA,CT3.pw.q.145=NA,CT4.pw.q.145=NA,CT5.pw.q.145=NA,CT6.pw.q.145=NA,CT1.pw.q.144=NA,CT2.pw.q.144=NA,CT3.pw.q.144=NA,CT4.pw.q.144=NA,CT5.pw.q.144=NA,CT6.pw.q.144=NA,CT1.pw.q.143=NA,CT2.pw.q.143=NA,
CT3.pw.q.143=NA,CT4.pw.q.143=NA,CT5.pw.q.143=NA,CT6.pw.q.143=NA,CT1.pw.q.142=NA,CT2.pw.q.142=NA,CT3.pw.q.142=NA,CT4.pw.q.142=NA,CT5.pw.q.142=NA,CT6.pw.q.142=NA,CT1.pw.q.141=NA,CT2.pw.q.141=NA,CT3.pw.q.141=NA,CT4.pw.q.141=NA,CT5.pw.q.141=NA,
CT6.pw.q.141=NA,CT1.pw.q.140=NA,CT2.pw.q.140=NA,CT3.pw.q.140=NA,CT4.pw.q.140=NA,CT5.pw.q.140=NA,CT6.pw.q.140=NA,CT1.pw.q.139=NA,CT2.pw.q.139=NA,CT3.pw.q.139=NA,CT4.pw.q.139=NA,CT5.pw.q.139=NA,CT6.pw.q.139=NA,CT1.pw.q.138=NA,CT2.pw.q.138=NA,
CT3.pw.q.138=NA,CT4.pw.q.138=NA,CT5.pw.q.138=NA,CT6.pw.q.138=NA,CT1.pw.q.137=NA,CT2.pw.q.137=NA,CT3.pw.q.137=NA,CT4.pw.q.137=NA,CT5.pw.q.137=NA,CT6.pw.q.137=NA,CT1.pw.q.136=NA,CT2.pw.q.136=NA,CT3.pw.q.136=NA,CT4.pw.q.136=NA,CT5.pw.q.136=NA,
CT6.pw.q.136=NA,CT1.pw.q.135=NA,CT2.pw.q.135=NA,CT3.pw.q.135=NA,CT4.pw.q.135=NA,CT5.pw.q.135=NA,CT6.pw.q.135=NA,CT1.pw.q.134=NA,CT2.pw.q.134=NA,CT3.pw.q.134=NA,CT4.pw.q.134=NA,CT5.pw.q.134=NA,CT6.pw.q.134=NA,CT1.pw.q.133=NA,CT2.pw.q.133=NA,
CT3.pw.q.133=NA,CT4.pw.q.133=NA,CT5.pw.q.133=NA,CT6.pw.q.133=NA,CT1.pw.q.132=NA,CT2.pw.q.132=NA,CT3.pw.q.132=NA,CT4.pw.q.132=NA,CT5.pw.q.132=NA,CT6.pw.q.132=NA,CT1.pw.q.131=NA,CT2.pw.q.131=NA,CT3.pw.q.131=NA,CT4.pw.q.131=NA,CT5.pw.q.131=NA,
CT6.pw.q.131=NA,CT1.pw.q.130=NA,CT2.pw.q.130=NA,CT3.pw.q.130=NA,CT4.pw.q.130=NA,CT5.pw.q.130=NA,CT6.pw.q.130=NA,CT1.pw.q.129=NA,CT2.pw.q.129=NA,CT3.pw.q.129=NA,CT4.pw.q.129=NA,CT5.pw.q.129=NA,CT6.pw.q.129=NA,CT1.pw.q.128=NA,CT2.pw.q.128=NA,
CT3.pw.q.128=NA,CT4.pw.q.128=NA,CT5.pw.q.128=NA,CT6.pw.q.128=NA,CT1.pw.q.127=NA,CT2.pw.q.127=NA,CT3.pw.q.127=NA,CT4.pw.q.127=NA,CT5.pw.q.127=NA,CT6.pw.q.127=NA,CT1.pw.q.126=NA,CT2.pw.q.126=NA,CT3.pw.q.126=NA,CT4.pw.q.126=NA,CT5.pw.q.126=NA,
CT6.pw.q.126=NA,CT1.pw.q.125=NA,CT2.pw.q.125=NA,CT3.pw.q.125=NA,CT4.pw.q.125=NA,CT5.pw.q.125=NA,CT6.pw.q.125=NA,CT1.pw.q.124=NA,CT2.pw.q.124=NA,CT3.pw.q.124=NA,CT4.pw.q.124=NA,CT5.pw.q.124=NA,CT6.pw.q.124=NA,CT1.pw.q.123=NA,CT2.pw.q.123=NA,
CT3.pw.q.123=NA,CT4.pw.q.123=NA,CT5.pw.q.123=NA,CT6.pw.q.123=NA,CT1.pw.q.122=NA,CT2.pw.q.122=NA,CT3.pw.q.122=NA,CT4.pw.q.122=NA,CT5.pw.q.122=NA,CT6.pw.q.122=NA,CT1.pw.q.121=NA,CT2.pw.q.121=NA,CT3.pw.q.121=NA,CT4.pw.q.121=NA,CT5.pw.q.121=NA,
CT6.pw.q.121=NA,CT1.pw.q.120=NA,CT2.pw.q.120=NA,CT3.pw.q.120=NA,CT4.pw.q.120=NA,CT5.pw.q.120=NA,CT6.pw.q.120=NA,CT1.pw.q.119=NA,CT2.pw.q.119=NA,CT3.pw.q.119=NA,CT4.pw.q.119=NA,CT5.pw.q.119=NA,CT6.pw.q.119=NA,CT1.pw.q.118=NA,CT2.pw.q.118=NA,
CT3.pw.q.118=NA,CT4.pw.q.118=NA,CT5.pw.q.118=NA,CT6.pw.q.118=NA,CT1.pw.q.117=NA,CT2.pw.q.117=NA,CT3.pw.q.117=NA,CT4.pw.q.117=NA,CT5.pw.q.117=NA,CT6.pw.q.117=NA,CT1.pw.q.116=NA,CT2.pw.q.116=NA,CT3.pw.q.116=NA,CT4.pw.q.116=NA,CT5.pw.q.116=NA,
CT6.pw.q.116=NA,CT1.pw.q.115=NA,CT2.pw.q.115=NA,CT3.pw.q.115=NA,CT4.pw.q.115=NA,CT5.pw.q.115=NA,CT6.pw.q.115=NA,CT1.pw.q.114=NA,CT2.pw.q.114=NA,CT3.pw.q.114=NA,CT4.pw.q.114=NA,CT5.pw.q.114=NA,CT6.pw.q.114=NA,CT1.pw.q.113=NA,CT2.pw.q.113=NA,
CT3.pw.q.113=NA,CT4.pw.q.113=NA,CT5.pw.q.113=NA,CT6.pw.q.113=NA,CT1.pw.q.112=NA,CT2.pw.q.112=NA,CT3.pw.q.112=NA,CT4.pw.q.112=NA,CT5.pw.q.112=NA,CT6.pw.q.112=NA,CT1.pw.q.111=NA,CT2.pw.q.111=NA,CT3.pw.q.111=NA,CT4.pw.q.111=NA,CT5.pw.q.111=NA,
CT6.pw.q.111=NA,CT1.pw.q.110=NA,CT2.pw.q.110=NA,CT3.pw.q.110=NA,CT4.pw.q.110=NA,CT5.pw.q.110=NA,CT6.pw.q.110=NA,CT1.pw.q.109=NA,CT2.pw.q.109=NA,CT3.pw.q.109=NA,CT4.pw.q.109=NA,CT5.pw.q.109=NA,CT6.pw.q.109=NA,CT1.pw.q.108=NA,CT2.pw.q.108=NA,
CT3.pw.q.108=NA,CT4.pw.q.108=NA,CT5.pw.q.108=NA,CT6.pw.q.108=NA,CT1.pw.q.107=NA,CT2.pw.q.107=NA,CT3.pw.q.107=NA,CT4.pw.q.107=NA,CT5.pw.q.107=NA,CT6.pw.q.107=NA,CT1.pw.q.106=NA,CT2.pw.q.106=NA,CT3.pw.q.106=NA,CT4.pw.q.106=NA,CT5.pw.q.106=NA,
CT6.pw.q.106=NA,CT1.pw.q.105=NA,CT2.pw.q.105=NA,CT3.pw.q.105=NA,CT4.pw.q.105=NA,CT5.pw.q.105=NA,CT6.pw.q.105=NA,CT1.pw.q.104=NA,CT2.pw.q.104=NA,CT3.pw.q.104=NA,CT4.pw.q.104=NA,CT5.pw.q.104=NA,CT6.pw.q.104=NA,CT1.pw.q.103=NA,CT2.pw.q.103=NA,
CT3.pw.q.103=NA,CT4.pw.q.103=NA,CT5.pw.q.103=NA,CT6.pw.q.103=NA,CT1.pw.q.102=NA,CT2.pw.q.102=NA,CT3.pw.q.102=NA,CT4.pw.q.102=NA,CT5.pw.q.102=NA,CT6.pw.q.102=NA,CT1.pw.q.101=NA,CT2.pw.q.101=NA,CT3.pw.q.101=NA,CT4.pw.q.101=NA,CT5.pw.q.101=NA,
CT6.pw.q.101=NA,CT1.pw.q.100=NA,CT2.pw.q.100=NA,CT3.pw.q.100=NA,CT4.pw.q.100=NA,CT5.pw.q.100=NA,CT6.pw.q.100=NA,CT1.pw.q.099=NA,CT2.pw.q.099=NA,CT3.pw.q.099=NA,CT4.pw.q.099=NA,CT5.pw.q.099=NA,CT6.pw.q.099=NA,CT1.pw.q.098=NA,CT2.pw.q.098=NA,
CT3.pw.q.098=NA,CT4.pw.q.098=NA,CT5.pw.q.098=NA,CT6.pw.q.098=NA,CT1.pw.q.097=NA,CT2.pw.q.097=NA,CT3.pw.q.097=NA,CT4.pw.q.097=NA,CT5.pw.q.097=NA,CT6.pw.q.097=NA,CT1.pw.q.096=NA,CT2.pw.q.096=NA,CT3.pw.q.096=NA,CT4.pw.q.096=NA,CT5.pw.q.096=NA,
CT6.pw.q.096=NA,CT1.pw.q.095=NA,CT2.pw.q.095=NA,CT3.pw.q.095=NA,CT4.pw.q.095=NA,CT5.pw.q.095=NA,CT6.pw.q.095=NA,CT1.pw.q.094=NA,CT2.pw.q.094=NA,CT3.pw.q.094=NA,CT4.pw.q.094=NA,CT5.pw.q.094=NA,CT6.pw.q.094=NA,CT1.pw.q.093=NA,CT2.pw.q.093=NA,
CT3.pw.q.093=NA,CT4.pw.q.093=NA,CT5.pw.q.093=NA,CT6.pw.q.093=NA,CT1.pw.q.092=NA,CT2.pw.q.092=NA,CT3.pw.q.092=NA,CT4.pw.q.092=NA,CT5.pw.q.092=NA,CT6.pw.q.092=NA,CT1.pw.q.091=NA,CT2.pw.q.091=NA,CT3.pw.q.091=NA,CT4.pw.q.091=NA,CT5.pw.q.091=NA,
CT6.pw.q.091=NA,CT1.pw.q.090=NA,CT2.pw.q.090=NA,CT3.pw.q.090=NA,CT4.pw.q.090=NA,CT5.pw.q.090=NA,CT6.pw.q.090=NA,CT1.pw.q.089=NA,CT2.pw.q.089=NA,CT3.pw.q.089=NA,CT4.pw.q.089=NA,CT5.pw.q.089=NA,CT6.pw.q.089=NA,CT1.pw.q.088=NA,CT2.pw.q.088=NA,
CT3.pw.q.088=NA,CT4.pw.q.088=NA,CT5.pw.q.088=NA,CT6.pw.q.088=NA,CT1.pw.q.087=NA,CT2.pw.q.087=NA,CT3.pw.q.087=NA,CT4.pw.q.087=NA,CT5.pw.q.087=NA,CT6.pw.q.087=NA,CT1.pw.q.086=NA,CT2.pw.q.086=NA,CT3.pw.q.086=NA,CT4.pw.q.086=NA,CT5.pw.q.086=NA,
CT6.pw.q.086=NA,CT1.pw.q.085=NA,CT2.pw.q.085=NA,CT3.pw.q.085=NA,CT4.pw.q.085=NA,CT5.pw.q.085=NA,CT6.pw.q.085=NA,CT1.pw.q.084=NA,CT2.pw.q.084=NA,CT3.pw.q.084=NA,CT4.pw.q.084=NA,CT5.pw.q.084=NA,CT6.pw.q.084=NA,CT1.pw.q.083=NA,CT2.pw.q.083=NA,
CT3.pw.q.083=NA,CT4.pw.q.083=NA,CT5.pw.q.083=NA,CT6.pw.q.083=NA,CT1.pw.q.082=NA,CT2.pw.q.082=NA,CT3.pw.q.082=NA,CT4.pw.q.082=NA,CT5.pw.q.082=NA,CT6.pw.q.082=NA,CT1.pw.q.081=NA,CT2.pw.q.081=NA,CT3.pw.q.081=NA,CT4.pw.q.081=NA,CT5.pw.q.081=NA,
CT6.pw.q.081=NA,CT1.pw.q.080=NA,CT2.pw.q.080=NA,CT3.pw.q.080=NA,CT4.pw.q.080=NA,CT5.pw.q.080=NA,CT6.pw.q.080=NA,CT1.pw.q.079=NA,CT2.pw.q.079=NA,CT3.pw.q.079=NA,CT4.pw.q.079=NA,CT5.pw.q.079=NA,CT6.pw.q.079=NA,CT1.pw.q.078=NA,CT2.pw.q.078=NA,
CT3.pw.q.078=NA,CT4.pw.q.078=NA,CT5.pw.q.078=NA,CT6.pw.q.078=NA,CT1.pw.q.077=NA,CT2.pw.q.077=NA,CT3.pw.q.077=NA,CT4.pw.q.077=NA,CT5.pw.q.077=NA,CT6.pw.q.077=NA,CT1.pw.q.076=NA,CT2.pw.q.076=NA,CT3.pw.q.076=NA,CT4.pw.q.076=NA,CT5.pw.q.076=NA,
CT6.pw.q.076=NA,CT1.pw.q.075=NA,CT2.pw.q.075=NA,CT3.pw.q.075=NA,CT4.pw.q.075=NA,CT5.pw.q.075=NA,CT6.pw.q.075=NA,CT1.pw.q.074=NA,CT2.pw.q.074=NA,CT3.pw.q.074=NA,CT4.pw.q.074=NA,CT5.pw.q.074=NA,CT6.pw.q.074=NA,CT1.pw.q.073=NA,CT2.pw.q.073=NA,
CT3.pw.q.073=NA,CT4.pw.q.073=NA,CT5.pw.q.073=NA,CT6.pw.q.073=NA,CT1.pw.q.072=NA,CT2.pw.q.072=NA,CT3.pw.q.072=NA,CT4.pw.q.072=NA,CT5.pw.q.072=NA,CT6.pw.q.072=NA,CT1.pw.q.071=NA,CT2.pw.q.071=NA,CT3.pw.q.071=NA,CT4.pw.q.071=NA,CT5.pw.q.071=NA,
CT6.pw.q.071=NA,CT1.pw.q.070=NA,CT2.pw.q.070=NA,CT3.pw.q.070=NA,CT4.pw.q.070=NA,CT5.pw.q.070=NA,CT6.pw.q.070=NA,CT1.pw.q.069=NA,CT2.pw.q.069=NA,CT3.pw.q.069=NA,CT4.pw.q.069=NA,CT5.pw.q.069=NA,CT6.pw.q.069=NA,CT1.pw.q.068=NA,CT2.pw.q.068=NA,
CT3.pw.q.068=NA,CT4.pw.q.068=NA,CT5.pw.q.068=NA,CT6.pw.q.068=NA,CT1.pw.q.067=NA,CT2.pw.q.067=NA,CT3.pw.q.067=NA,CT4.pw.q.067=NA,CT5.pw.q.067=NA,CT6.pw.q.067=NA,CT1.pw.q.066=NA,CT2.pw.q.066=NA,CT3.pw.q.066=NA,CT4.pw.q.066=NA,CT5.pw.q.066=NA,
CT6.pw.q.066=NA,CT1.pw.q.065=NA,CT2.pw.q.065=NA,CT3.pw.q.065=NA,CT4.pw.q.065=NA,CT5.pw.q.065=NA,CT6.pw.q.065=NA,CT1.pw.q.064=NA,CT2.pw.q.064=NA,CT3.pw.q.064=NA,CT4.pw.q.064=NA,CT5.pw.q.064=NA,CT6.pw.q.064=NA,CT1.pw.q.063=NA,CT2.pw.q.063=NA,
CT3.pw.q.063=NA,CT4.pw.q.063=NA,CT5.pw.q.063=NA,CT6.pw.q.063=NA,CT1.pw.q.062=NA,CT2.pw.q.062=NA,CT3.pw.q.062=NA,CT4.pw.q.062=NA,CT5.pw.q.062=NA,CT6.pw.q.062=NA,CT1.pw.q.061=NA,CT2.pw.q.061=NA,CT3.pw.q.061=NA,CT4.pw.q.061=NA,CT5.pw.q.061=NA,
CT6.pw.q.061=NA,CT1.pw.q.060=NA,CT2.pw.q.060=NA,CT3.pw.q.060=NA,CT4.pw.q.060=NA,CT5.pw.q.060=NA,CT6.pw.q.060=NA,CT1.pw.q.059=NA,CT2.pw.q.059=NA,CT3.pw.q.059=NA,CT4.pw.q.059=NA,CT5.pw.q.059=NA,CT6.pw.q.059=NA,CT1.pw.q.058=NA,CT2.pw.q.058=NA,
CT3.pw.q.058=NA,CT4.pw.q.058=NA,CT5.pw.q.058=NA,CT6.pw.q.058=NA,CT1.pw.q.057=NA,CT2.pw.q.057=NA,CT3.pw.q.057=NA,CT4.pw.q.057=NA,CT5.pw.q.057=NA,CT6.pw.q.057=NA,CT1.pw.q.056=NA,CT2.pw.q.056=NA,CT3.pw.q.056=NA,CT4.pw.q.056=NA,CT5.pw.q.056=NA,
CT6.pw.q.056=NA,CT1.pw.q.055=NA,CT2.pw.q.055=NA,CT3.pw.q.055=NA,CT4.pw.q.055=NA,CT5.pw.q.055=NA,CT6.pw.q.055=NA,CT1.pw.q.054=NA,CT2.pw.q.054=NA,CT3.pw.q.054=NA,CT4.pw.q.054=NA,CT5.pw.q.054=NA,CT6.pw.q.054=NA,CT1.pw.q.053=NA,CT2.pw.q.053=NA,
CT3.pw.q.053=NA,CT4.pw.q.053=NA,CT5.pw.q.053=NA,CT6.pw.q.053=NA,CT1.pw.q.052=NA,CT2.pw.q.052=NA,CT3.pw.q.052=NA,CT4.pw.q.052=NA,CT5.pw.q.052=NA,CT6.pw.q.052=NA,CT1.pw.q.051=NA,CT2.pw.q.051=NA,CT3.pw.q.051=NA,CT4.pw.q.051=NA,CT5.pw.q.051=NA,
CT6.pw.q.051=NA,CT1.pw.q.050=NA,CT2.pw.q.050=NA,CT3.pw.q.050=NA,CT4.pw.q.050=NA,CT5.pw.q.050=NA,CT6.pw.q.050=NA,CT1.pw.q.049=NA,CT2.pw.q.049=NA,CT3.pw.q.049=NA,CT4.pw.q.049=NA,CT5.pw.q.049=NA,CT6.pw.q.049=NA,CT1.pw.q.048=NA,CT2.pw.q.048=NA,
CT3.pw.q.048=NA,CT4.pw.q.048=NA,CT5.pw.q.048=NA,CT6.pw.q.048=NA,CT1.pw.q.047=NA,CT2.pw.q.047=NA,CT3.pw.q.047=NA,CT4.pw.q.047=NA,CT5.pw.q.047=NA,CT6.pw.q.047=NA,CT1.pw.q.046=NA,CT2.pw.q.046=NA,CT3.pw.q.046=NA,CT4.pw.q.046=NA,CT5.pw.q.046=NA,
CT6.pw.q.046=NA,CT1.pw.q.045=NA,CT2.pw.q.045=NA,CT3.pw.q.045=NA,CT4.pw.q.045=NA,CT5.pw.q.045=NA,CT6.pw.q.045=NA,CT1.pw.q.044=NA,CT2.pw.q.044=NA,CT3.pw.q.044=NA,CT4.pw.q.044=NA,CT5.pw.q.044=NA,CT6.pw.q.044=NA,CT1.pw.q.043=NA,CT2.pw.q.043=NA,
CT3.pw.q.043=NA,CT4.pw.q.043=NA,CT5.pw.q.043=NA,CT6.pw.q.043=NA,CT1.pw.q.042=NA,CT2.pw.q.042=NA,CT3.pw.q.042=NA,CT4.pw.q.042=NA,CT5.pw.q.042=NA,CT6.pw.q.042=NA,CT1.pw.q.041=NA,CT2.pw.q.041=NA,CT3.pw.q.041=NA,CT4.pw.q.041=NA,CT5.pw.q.041=NA,
CT6.pw.q.041=NA,CT1.pw.q.040=NA,CT2.pw.q.040=NA,CT3.pw.q.040=NA,CT4.pw.q.040=NA,CT5.pw.q.040=NA,CT6.pw.q.040=NA,CT1.pw.q.039=NA,CT2.pw.q.039=NA,CT3.pw.q.039=NA,CT4.pw.q.039=NA,CT5.pw.q.039=NA,CT6.pw.q.039=NA,CT1.pw.q.038=NA,CT2.pw.q.038=NA,
CT3.pw.q.038=NA,CT4.pw.q.038=NA,CT5.pw.q.038=NA,CT6.pw.q.038=NA,CT1.pw.q.037=NA,CT2.pw.q.037=NA,CT3.pw.q.037=NA,CT4.pw.q.037=NA,CT5.pw.q.037=NA,CT6.pw.q.037=NA,CT1.pw.q.036=NA,CT2.pw.q.036=NA,CT3.pw.q.036=NA,CT4.pw.q.036=NA,CT5.pw.q.036=NA,
CT6.pw.q.036=NA,CT1.pw.q.035=NA,CT2.pw.q.035=NA,CT3.pw.q.035=NA,CT4.pw.q.035=NA,CT5.pw.q.035=NA,CT6.pw.q.035=NA,CT1.pw.q.034=NA,CT2.pw.q.034=NA,CT3.pw.q.034=NA,CT4.pw.q.034=NA,CT5.pw.q.034=NA,CT6.pw.q.034=NA,CT1.pw.q.033=NA,CT2.pw.q.033=NA,
CT3.pw.q.033=NA,CT4.pw.q.033=NA,CT5.pw.q.033=NA,CT6.pw.q.033=NA,CT1.pw.q.032=NA,CT2.pw.q.032=NA,CT3.pw.q.032=NA,CT4.pw.q.032=NA,CT5.pw.q.032=NA,CT6.pw.q.032=NA,CT1.pw.q.031=NA,CT2.pw.q.031=NA,CT3.pw.q.031=NA,CT4.pw.q.031=NA,CT5.pw.q.031=NA,
CT6.pw.q.031=NA,CT1.pw.q.030=NA,CT2.pw.q.030=NA,CT3.pw.q.030=NA,CT4.pw.q.030=NA,CT5.pw.q.030=NA,CT6.pw.q.030=NA,CT1.pw.q.029=NA,CT2.pw.q.029=NA,CT3.pw.q.029=NA,CT4.pw.q.029=NA,CT5.pw.q.029=NA,CT6.pw.q.029=NA,CT1.pw.q.028=NA,CT2.pw.q.028=NA,
CT3.pw.q.028=NA,CT4.pw.q.028=NA,CT5.pw.q.028=NA,CT6.pw.q.028=NA,CT1.pw.q.027=NA,CT2.pw.q.027=NA,CT3.pw.q.027=NA,CT4.pw.q.027=NA,CT5.pw.q.027=NA,CT6.pw.q.027=NA,CT1.pw.q.026=NA,CT2.pw.q.026=NA,CT3.pw.q.026=NA,CT4.pw.q.026=NA,CT5.pw.q.026=NA,
CT6.pw.q.026=NA,CT1.pw.q.025=NA,CT2.pw.q.025=NA,CT3.pw.q.025=NA,CT4.pw.q.025=NA,CT5.pw.q.025=NA,CT6.pw.q.025=NA,CT1.pw.q.024=NA,CT2.pw.q.024=NA,CT3.pw.q.024=NA,CT4.pw.q.024=NA,CT5.pw.q.024=NA,CT6.pw.q.024=NA,CT1.pw.q.023=NA,CT2.pw.q.023=NA,
CT3.pw.q.023=NA,CT4.pw.q.023=NA,CT5.pw.q.023=NA,CT6.pw.q.023=NA,CT1.pw.q.022=NA,CT2.pw.q.022=NA,CT3.pw.q.022=NA,CT4.pw.q.022=NA,CT5.pw.q.022=NA,CT6.pw.q.022=NA,CT1.pw.q.021=NA,CT2.pw.q.021=NA,CT3.pw.q.021=NA,CT4.pw.q.021=NA,CT5.pw.q.021=NA,
CT6.pw.q.021=NA,CT1.pw.q.020=NA,CT2.pw.q.020=NA,CT3.pw.q.020=NA,CT4.pw.q.020=NA,CT5.pw.q.020=NA,CT6.pw.q.020=NA,CT1.pw.q.019=NA,CT2.pw.q.019=NA,CT3.pw.q.019=NA,CT4.pw.q.019=NA,CT5.pw.q.019=NA,CT6.pw.q.019=NA,CT1.pw.q.018=NA,CT2.pw.q.018=NA,
CT3.pw.q.018=NA,CT4.pw.q.018=NA,CT5.pw.q.018=NA,CT6.pw.q.018=NA,CT1.pw.q.017=NA,CT2.pw.q.017=NA,CT3.pw.q.017=NA,CT4.pw.q.017=NA,CT5.pw.q.017=NA,CT6.pw.q.017=NA,CT1.pw.q.016=NA,CT2.pw.q.016=NA,CT3.pw.q.016=NA,CT4.pw.q.016=NA,CT5.pw.q.016=NA,
CT6.pw.q.016=NA,CT1.pw.q.015=NA,CT2.pw.q.015=NA,CT3.pw.q.015=NA,CT4.pw.q.015=NA,CT5.pw.q.015=NA,CT6.pw.q.015=NA,CT1.pw.q.014=NA,CT2.pw.q.014=NA,CT3.pw.q.014=NA,CT4.pw.q.014=NA,CT5.pw.q.014=NA,CT6.pw.q.014=NA,CT1.pw.q.013=NA,CT2.pw.q.013=NA,
CT3.pw.q.013=NA,CT4.pw.q.013=NA,CT5.pw.q.013=NA,CT6.pw.q.013=NA,CT1.pw.q.012=NA,CT2.pw.q.012=NA,CT3.pw.q.012=NA,CT4.pw.q.012=NA,CT5.pw.q.012=NA,CT6.pw.q.012=NA,CT1.pw.q.011=NA,CT2.pw.q.011=NA,CT3.pw.q.011=NA,CT4.pw.q.011=NA,CT5.pw.q.011=NA,
CT6.pw.q.011=NA,CT1.pw.q.010=NA,CT2.pw.q.010=NA,CT3.pw.q.010=NA,CT4.pw.q.010=NA,CT5.pw.q.010=NA,CT6.pw.q.010=NA,CT1.pw.q.009=NA,CT2.pw.q.009=NA,CT3.pw.q.009=NA,CT4.pw.q.009=NA,CT5.pw.q.009=NA,CT6.pw.q.009=NA,CT1.pw.q.008=NA,CT2.pw.q.008=NA,
CT3.pw.q.008=NA,CT4.pw.q.008=NA,CT5.pw.q.008=NA,CT6.pw.q.008=NA,CT1.pw.q.007=NA,CT2.pw.q.007=NA,CT3.pw.q.007=NA,CT4.pw.q.007=NA,CT5.pw.q.007=NA,CT6.pw.q.007=NA,CT1.pw.q.006=NA,CT2.pw.q.006=NA,CT3.pw.q.006=NA,CT4.pw.q.006=NA,CT5.pw.q.006=NA,
CT6.pw.q.006=NA,CT1.pw.q.005=NA,CT2.pw.q.005=NA,CT3.pw.q.005=NA,CT4.pw.q.005=NA,CT5.pw.q.005=NA,CT6.pw.q.005=NA,CT1.pw.q.004=NA,CT2.pw.q.004=NA,CT3.pw.q.004=NA,CT4.pw.q.004=NA,CT5.pw.q.004=NA,CT6.pw.q.004=NA,CT1.pw.q.003=NA,CT2.pw.q.003=NA,
CT3.pw.q.003=NA,CT4.pw.q.003=NA,CT5.pw.q.003=NA,CT6.pw.q.003=NA,CT1.pw.q.002=NA,CT2.pw.q.002=NA,CT3.pw.q.002=NA,CT4.pw.q.002=NA,CT5.pw.q.002=NA,CT6.pw.q.002=NA,CT1.pw.q.001=NA,CT2.pw.q.001=NA,CT3.pw.q.001=NA,CT4.pw.q.001=NA,CT5.pw.q.001=NA,
CT6.pw.q.001=NA
)
results = results[-1,]

for(rrr in 1:length(supercl)){ ######################################################################################################

ctime = Sys.time()

deltas = sample( c(1,1,1,1,1,1) ) * delta.pow
nfacts = sample( c(-1,-1,-1,1,1,1) )
nidxs = nfacts < 0

concentration_parameters = rev(c(0.2,0.4,0.8,1.6,3.2,6.4)) 

cell_fractions = rdirichlet(n=num_subjects,alpha=concentration_parameters)
cft = t(cell_fractions)
num_ctypes = ncol(cell_fractions)

cell_fractions_q = rdirichlet(n=num_subjects,alpha=concentration_parameters)
cft_q = t(cell_fractions_q)

clusters = supercl[[rrr]]
max_clst_length = max(unlist(lapply(clusters,FUN=length)))

cluster_info_tab = data.frame(pos=NA,cluster=NA)
for(clst in 1:length(clusters)){
	cluster = clusters[[clst]]
	for(cpg in 1:length(cluster) ){
		cluster_info_tab = rbind(cluster_info_tab,data.frame(pos=cluster[cpg],cluster=clst))
	}
}
cluster_info_tab = cluster_info_tab[-1,]

delta_mat = matrix(NA,nrow=6,ncol=nrow(cluster_info_tab))
delta_cor_mat = matrix(NA,nrow=nrow(cluster_info_tab),ncol=nrow(cluster_info_tab))
for(idx1 in 1:nrow(cluster_info_tab)){
	for(idx2 in 1:nrow(cluster_info_tab)){
		if(idx1 == idx2){
			delta_cor_mat[idx1,idx2] = 1
		} else if(cluster_info_tab[idx1,"cluster"]==cluster_info_tab[idx2,"cluster"]){
			delta_cor_mat[idx1,idx2] = cor.c.w
		} else{
			delta_cor_mat[idx1,idx2] = cor.c.b
		}
	}
}
for(idx3 in 1:6){
	delta_mat[idx3,] = mvrnorm(n=1, mu=rep(deltas[idx3],nrow(cluster_info_tab)), Sigma=var.cluster*delta_cor_mat)
}
delta_list = list()
counter = 1
for(clst in 1:length(clusters)){
	delta_list[[clst]] = list()
	for(idx4 in 1:length(clusters[[clst]])){
		delta_list[[clst]][[idx4]] = delta_mat[,counter]
		counter = counter+1
	}
}

pheno.sub = pheno.normal6[draw_random_mixture_candidates(pheno.normal6),]
selected_samples = as.character(pheno.sub$Array)
betas_sub = betas.normal6[,selected_samples ]

GRPS = rep(0:1,each=num_subjects/2)
X = array(NA,dim=c(num_subjects,length(clusters),max_clst_length))
QPRED = array(NA,dim=c(num_subjects,length(clusters),max_clst_length))
num_cpgs = c()


cat(
	"processing super cluster",rrr,"/",length(supercl),"(",paste(sort(unlist(lapply(clusters,FUN=length))),collapse=","),") with",length(unlist(clusters)),
	"cpgs ; delta.pow =",delta.pow, "; sd(x,z,c) =",paste0("(",noise.x,",",noise.z,",",noise.c,")"),"; cor(w,b) =",paste0("(",cor.c.w,",",cor.c.b,")"),
	"...\n"
)
	
	allAN = anchr[anchr$pos %in% unlist(clusters),]
	allBE = betas_sub[row.names(betas_sub) %in% row.names(allAN),]
	ct_median = c()
	ct_mins = c()
	ct_maxs = c()
	cases = c()
	for(ppp in 1:num_ctypes){
		ct_median[ppp] = median(allBE[,ppp])
		ct_mins[ppp] = min(allBE[,ppp])
		ct_maxs[ppp] = max(allBE[,ppp])
		
		if(ct_maxs[ppp] < 0.6){
			cases[ppp] = 1
		} else if(ct_mins[ppp] > 0.4){
			cases[ppp] = 2
		} else{
			if(ct_median[ppp] < 0.5){
				cases[ppp] = 1
			} else{ cases[ppp] = 2 }
		}
	}
	
for(clst in 1:length(clusters)){
	TNCPG = length(clusters[[clst]])
	subAN = anchr[anchr$pos %in% clusters[[clst]],]
	subBE = betas_sub[row.names(betas_sub) %in% row.names(subAN),]
	subMixBE = matrix(NA,nrow=TNCPG,ncol=num_subjects)
	
	for(si in 1:num_subjects){
		tmpBE = matrix(subBE,nrow=TNCPG)
		for(sss1 in 1:TNCPG){ # CPG INDEX
			for(sss2 in 1:num_ctypes){ # CELL TYPE INDEX
				delta.gam = c(0,delta_list[[clst]][[sss1]][sss2])
				if(cases[sss2] == 2){
					delta.gam = c(-delta_list[[clst]][[sss1]][sss2],0)
				}
				tmpBE[sss1,sss2] = tmpBE[sss1,sss2] + delta.gam[GRPS[si]+1]*nfacts[sss2] + rnorm(n=1,mean=0,sd=noise.z)
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
				delta.gam = c(0,delta_list[[clst]][[sss1]][sss2])
				if(cases[sss2] == 2){
					delta.gam = c(-delta_list[[clst]][[sss1]][sss2],0)
				}
				tmpBE[sss1,sss2] = tmpBE[sss1,sss2] + delta.gam[GRPS[si]+1]*nfacts[sss2] + rnorm(n=1,mean=0,sd=noise.z)
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

inivals = generate.initial.values.berry.m1(clusters=clusters, chains=1, KCF=6, N=num_subjects)

ftime = Sys.time()
fit_berry_m2 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","wrapper.gam","tau","psi","omu","oomu","ogam","oogam",
		"otau","ootau","olambda","oolambda", 
		"RMSE","PREDRMSE"
	), 
	model.file = berry_like_model_m2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)
fit_times[["Berry.M2"]] = as.numeric(difftime(Sys.time(),ftime,units = "secs"))

ftime = Sys.time()
fit_berry_m2_cprior <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","wrapper.gam","tau","psi","omu","oomu","ogam","oogam",
		"otau","ootau","olambda","oolambda", 
		"RMSE","PREDRMSE"
	), 
	model.file = berry_like_model_m2_custom_prior, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)
fit_times[["Berry.M2.cprior1"]] = as.numeric(difftime(Sys.time(),ftime,units = "secs"))

ftime = Sys.time()
fit_berry_m2_cprior2 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"wrapper.mu","wrapper.gam","tau","psi","omu","oomu","ogam","oogam",
		"otau","ootau","olambda","oolambda", 
		"RMSE","PREDRMSE"
	), 
	model.file = berry_like_model_m2_custom_prior2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)
fit_times[["Berry.M2.cprior2"]] = as.numeric(difftime(Sys.time(),ftime,units = "secs"))

ftime = Sys.time()
fit_tca <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"mu","gam","tau","psi","RMSE","PREDRMSE"
	), 
	model.file = tca_like_model_flip, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)
fit_times[["TCA"]] = as.numeric(difftime(Sys.time(),ftime,units = "secs"))

ftime = Sys.time()
fit_tca_cprior <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"mu","gam","tau","psi","RMSE","PREDRMSE"
	), 
	model.file = tca_like_model_flip_custom_prior, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)
fit_times[["TCA.cprior1"]] = as.numeric(difftime(Sys.time(),ftime,units = "secs"))

ftime = Sys.time()
fit_tca_cprior2 <- bugs(
	dats, inits = inivals, 
	parameters.to.save = c(
		"mu","gam","tau","psi","RMSE","PREDRMSE"
	), 
	model.file = tca_like_model_flip_custom_prior2, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8
)
fit_times[["TCA.cprior2"]] = as.numeric(difftime(Sys.time(),ftime,units = "secs"))


qrejects0 = c()
qrejects0cp = c()
qrejects0cp2 = c()
qrejectsBerryM2 = c()
qrejectsBerryM2cp = c()
qrejectsBerryM2cp2 = c()

indi_qrejects0 = matrix(NA,ncol=1, nrow=1)
indi_qrejects0cp = matrix(NA,ncol=1, nrow=1)
indi_qrejects0cp2 = matrix(NA,ncol=1, nrow=1)
indi_qrejectsBerryM2 = matrix(NA,ncol=1, nrow=1)
indi_qrejectsBerryM2cp = matrix(NA,ncol=1, nrow=1)
indi_qrejectsBerryM2cp2 = matrix(NA,ncol=1, nrow=1)

for(zz in 1:length(tquants)){
	qqs0 = c()
	qqs0cp = c()
	qqs0cp2 = c()
	qqs_berry_m2 = c()
	qqs_berry_m2cp = c()
	qqs_berry_m2cp2 = c()
	
	qqs0_indi = matrix(NA,ncol=6, nrow=0)
	qqs0cp_indi = matrix(NA,ncol=6, nrow=0)
	qqs0cp2_indi = matrix(NA,ncol=6, nrow=0)
	qqs_berry_m2_indi = matrix(NA,ncol=6, nrow=0)
	qqs_berry_m2cp_indi = matrix(NA,ncol=6, nrow=0)
	qqs_berry_m2cp2_indi = matrix(NA,ncol=6, nrow=0)
	
	for(zc in 1:length(clusters)){
		for(zi in 1:num_ctypes){
		
			ctquant = tquants[zz]
			
			for(zj in 1:length(clusters[[zc]])){
				qqs0 = c( qqs0, (quantile(fit_tca$sims.list[["gam"]][,zi,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_tca$sims.list[["gam"]][,zi,zc,zj], probs=1-ctquant/2) < 0) )
				qqs0cp = c( qqs0cp, (quantile(fit_tca_cprior$sims.list[["gam"]][,zi,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_tca_cprior$sims.list[["gam"]][,zi,zc,zj], probs=1-ctquant/2) < 0) )
				qqs0cp2 = c( qqs0cp2, (quantile(fit_tca_cprior2$sims.list[["gam"]][,zi,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_tca_cprior2$sims.list[["gam"]][,zi,zc,zj], probs=1-ctquant/2) < 0) )
				qqs_berry_m2 = c( qqs_berry_m2, (quantile(fit_berry_m2$sims.list[["wrapper.gam"]][,zi,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_berry_m2$sims.list[["wrapper.gam"]][,zi,zc,zj], probs=1-ctquant/2) < 0) )
				qqs_berry_m2cp = c( qqs_berry_m2cp, (quantile(fit_berry_m2_cprior$sims.list[["wrapper.gam"]][,zi,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_berry_m2_cprior$sims.list[["wrapper.gam"]][,zi,zc,zj], probs=1-ctquant/2) < 0) )
				qqs_berry_m2cp2 = c( qqs_berry_m2cp2, (quantile(fit_berry_m2_cprior2$sims.list[["wrapper.gam"]][,zi,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_berry_m2_cprior2$sims.list[["wrapper.gam"]][,zi,zc,zj], probs=1-ctquant/2) < 0) )
			}
		}
		
		for(zj in 1:length(clusters[[zc]])){
			qqs0_itrow = c()
			qqs0cp_itrow = c()
			qqs0cp2_itrow = c()
			qqs_berry_m2_itrow = c()
			qqs_berry_m2cp_itrow = c()
			qqs_berry_m2cp2_itrow = c()
			for(zcur in 1:6){
			
				ctquant = tquants[zz]
				
				qqs0_itrow[zcur] = (quantile(fit_tca$sims.list[["gam"]][,zcur,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_tca$sims.list[["gam"]][,zcur,zc,zj], probs=1-ctquant/2) < 0)
				qqs0cp_itrow[zcur] = (quantile(fit_tca_cprior$sims.list[["gam"]][,zcur,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_tca_cprior$sims.list[["gam"]][,zcur,zc,zj], probs=1-ctquant/2) < 0)
				qqs0cp2_itrow[zcur] = (quantile(fit_tca_cprior2$sims.list[["gam"]][,zcur,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_tca_cprior2$sims.list[["gam"]][,zcur,zc,zj], probs=1-ctquant/2) < 0)
				qqs_berry_m2_itrow[zcur] = (quantile(fit_berry_m2$sims.list[["wrapper.gam"]][,zcur,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_berry_m2$sims.list[["wrapper.gam"]][,zcur,zc,zj], probs=1-ctquant/2) < 0)
				qqs_berry_m2cp_itrow[zcur] = (quantile(fit_berry_m2_cprior$sims.list[["wrapper.gam"]][,zcur,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_berry_m2_cprior$sims.list[["wrapper.gam"]][,zcur,zc,zj], probs=1-ctquant/2) < 0)
				qqs_berry_m2cp2_itrow[zcur] = (quantile(fit_berry_m2_cprior2$sims.list[["wrapper.gam"]][,zcur,zc,zj], probs=ctquant/2) > 0) || (quantile(fit_berry_m2_cprior2$sims.list[["wrapper.gam"]][,zcur,zc,zj], probs=1-ctquant/2) < 0)
			}
			qqs0_indi = rbind(qqs0_indi, qqs0_itrow)
			qqs0cp_indi = rbind(qqs0cp_indi, qqs0cp_itrow)
			qqs0cp2_indi = rbind(qqs0cp2_indi, qqs0cp2_itrow)
			qqs_berry_m2_indi = rbind(qqs_berry_m2_indi, qqs_berry_m2_itrow)
			qqs_berry_m2cp_indi = rbind(qqs_berry_m2cp_indi, qqs_berry_m2cp_itrow)
			qqs_berry_m2cp2_indi = rbind(qqs_berry_m2cp2_indi, qqs_berry_m2cp2_itrow)
		}
	}
	qrejects0[zz] = mean(qqs0>0)
	qrejects0cp[zz] = mean(qqs0cp>0)
	qrejects0cp2[zz] = mean(qqs0cp2>0)
	qrejectsBerryM2[zz] = mean(qqs_berry_m2>0)
	qrejectsBerryM2cp[zz] = mean(qqs_berry_m2cp>0)
	qrejectsBerryM2cp2[zz] = mean(qqs_berry_m2cp2>0)
	
	indi_qrejects0 = cbind( indi_qrejects0, mean(qqs0_indi[,1]>0),mean(qqs0_indi[,2]>0),mean(qqs0_indi[,3]>0),mean(qqs0_indi[,4]>0),mean(qqs0_indi[,5]>0),mean(qqs0_indi[,6]>0) )
	indi_qrejects0cp = cbind( indi_qrejects0cp, mean(qqs0cp_indi[,1]>0),mean(qqs0cp_indi[,2]>0),mean(qqs0cp_indi[,3]>0),mean(qqs0cp_indi[,4]>0),mean(qqs0cp_indi[,5]>0),mean(qqs0cp_indi[,6]>0) )
	indi_qrejects0cp2 = cbind( indi_qrejects0cp2, mean(qqs0cp2_indi[,1]>0),mean(qqs0cp2_indi[,2]>0),mean(qqs0cp2_indi[,3]>0),mean(qqs0cp2_indi[,4]>0),mean(qqs0cp2_indi[,5]>0),mean(qqs0cp2_indi[,6]>0) )
	indi_qrejectsBerryM2 = cbind( indi_qrejectsBerryM2, mean(qqs_berry_m2_indi[,1]>0),mean(qqs_berry_m2_indi[,2]>0),mean(qqs_berry_m2_indi[,3]>0),mean(qqs_berry_m2_indi[,4]>0),mean(qqs_berry_m2_indi[,5]>0),mean(qqs_berry_m2_indi[,6]>0) )
	indi_qrejectsBerryM2cp = cbind( indi_qrejectsBerryM2cp, mean(qqs_berry_m2cp_indi[,1]>0),mean(qqs_berry_m2cp_indi[,2]>0),mean(qqs_berry_m2cp_indi[,3]>0),mean(qqs_berry_m2cp_indi[,4]>0),mean(qqs_berry_m2cp_indi[,5]>0),mean(qqs_berry_m2cp_indi[,6]>0) )
	indi_qrejectsBerryM2cp2 = cbind( indi_qrejectsBerryM2cp2, mean(qqs_berry_m2cp2_indi[,1]>0),mean(qqs_berry_m2cp2_indi[,2]>0),mean(qqs_berry_m2cp2_indi[,3]>0),mean(qqs_berry_m2cp2_indi[,4]>0),mean(qqs_berry_m2cp2_indi[,5]>0),mean(qqs_berry_m2cp2_indi[,6]>0) )
}
indi_qrejects0 = indi_qrejects0[,-1]
indi_qrejects0cp = indi_qrejects0cp[,-1]
indi_qrejects0cp2 = indi_qrejects0cp2[,-1]
indi_qrejectsBerryM2 = indi_qrejectsBerryM2[,-1]
indi_qrejectsBerryM2cp = indi_qrejectsBerryM2cp[,-1]
indi_qrejectsBerryM2cp2 = indi_qrejectsBerryM2cp2[,-1]

QREJECTS = rbind(
	qrejects0, qrejects0cp, qrejects0cp2, 
	qrejectsBerryM2, qrejectsBerryM2cp, qrejectsBerryM2cp2
)
QREJECTS.CNT = QREJECTS * length(unlist(clusters)) * 6
colnames(QREJECTS) = paste0("pw.q.",sprintf("%03.f", tquants*1000))
colnames(QREJECTS.CNT) = paste0("cnt.q.",sprintf("%03.f", tquants*1000))

QindiREJECTS = rbind(
	indi_qrejects0, indi_qrejects0cp, indi_qrejects0cp2,
	indi_qrejectsBerryM2, indi_qrejectsBerryM2cp, indi_qrejectsBerryM2cp2
)
colnames(QindiREJECTS) = paste0(
"CT",
rep(1:6,times=length(tquants)),
".",
rep(
	paste0("pw.q.",sprintf("%03.f", tquants*1000)), each=6
)
)

i.res.tab = data.frame(
	IT.ID = rrr, MCMC.NUM=MCMC_NUM, MCMC.BURN=MCMC_BURN, MAX.SCN=MAX_SCN, MAX.CN=MAX_CN,
	noise.z=noise.z, noise.x=noise.x, noise.c=noise.c, cor.c.w=cor.c.w, cor.c.b = cor.c.b,
	model = c(
		"TCA", "TCA.cprior1", "TCA.cprior2", 
		"Berry.M2", "Berry.M2.cprior1", "Berry.M2.cprior2"
	),
	fit.time = c(
		fit_times[["TCA"]], fit_times[["TCA.cprior1"]], fit_times[["TCA.cprior2"]], 
		fit_times[["Berry.M2"]], fit_times[["Berry.M2.cprior1"]], fit_times[["Berry.M2.cprior2"]]
	),
	num.clusters = length(clusters),
	num.cpgs = length(unlist(clusters)),
	cluster.sizes = paste(sort(unlist(lapply(clusters,FUN=length))),collapse=","),
	super.cluster.span = calculate_super_cluster_span(clusters),
	ncr.dev = c(
		cross.rate(fit_tca$sims.list[["deviance"]]),
		cross.rate(fit_tca_cprior$sims.list[["deviance"]]),
		cross.rate(fit_tca_cprior2$sims.list[["deviance"]]),
		cross.rate(fit_berry_m2$sims.list[["deviance"]]),
		cross.rate(fit_berry_m2_cprior$sims.list[["deviance"]]),
		cross.rate(fit_berry_m2_cprior2$sims.list[["deviance"]])
		
		# deprecated autocorrelation code
		### mean(sapply(1:10,chain=fit_tca$sims.list[["deviance"]],FUN=autocor)),
	),
	PD = c(
		fit_tca$pD, fit_tca_cprior$pD, fit_tca_cprior2$pD, 
		fit_berry_m2$pD, fit_berry_m2_cprior$pD, fit_berry_m2_cprior2$pD
	),
	DIC = c(
		fit_tca$DIC, fit_tca_cprior$DIC, fit_tca_cprior2$DIC, 
		fit_berry_m2$DIC, fit_berry_m2_cprior$DIC, fit_berry_m2_cprior2$DIC
	),
	RMSE = c(
		fit_tca$summary[row.names(fit_tca$summary) == "RMSE",1], 
		fit_tca_cprior$summary[row.names(fit_tca_cprior$summary) == "RMSE",1], 
		fit_tca_cprior2$summary[row.names(fit_tca_cprior2$summary) == "RMSE",1], 
		fit_berry_m2$summary[row.names(fit_berry_m2$summary) == "RMSE",1], 
		fit_berry_m2_cprior$summary[row.names(fit_berry_m2_cprior$summary) == "RMSE",1], 
		fit_berry_m2_cprior2$summary[row.names(fit_berry_m2_cprior2$summary) == "RMSE",1]
	),
	PREDRMSE = c(
		fit_tca$summary[row.names(fit_tca$summary) == "PREDRMSE",1], 
		fit_tca_cprior$summary[row.names(fit_tca_cprior$summary) == "PREDRMSE",1], 
		fit_tca_cprior2$summary[row.names(fit_tca_cprior2$summary) == "PREDRMSE",1], 
		fit_berry_m2$summary[row.names(fit_berry_m2$summary) == "PREDRMSE",1], 
		fit_berry_m2_cprior$summary[row.names(fit_berry_m2_cprior$summary) == "PREDRMSE",1], 
		fit_berry_m2_cprior2$summary[row.names(fit_berry_m2_cprior2$summary) == "PREDRMSE",1]
	),
	QREJECTS, QREJECTS.CNT, QindiREJECTS
)

results = rbind(results,i.res.tab)

out_string = paste0(
	outpath,
	"two_arm_simC_cth1.",cluster_threshold1,"_cth2.",cluster_threshold2,"_",chrom,
	"_nsbj.",num_subjects,"_mcn.",MAX_CN,"_mscn.",MAX_SCN,"_delta.",delta.pow,
	"_noise.z.",noise.z,"_noise.x.",noise.x,"_noise.c.",noise.c,"_cor.c.w.",cor.c.w,"_cor.c.b.",cor.c.b,
	"_table.csv"
)
write.csv(results,file=out_string,row.names=FALSE)

cat("   ... finished iteration after a "); print(Sys.time()-ctime)

rm(
	dats, X, QPRED, cell_fractions, cell_fractions_q, cft, cft_q,
	fit_tca,fit_tca_cprior,fit_tca_cprior2,fit_tca_cprior3,
	fit_berry_m2,fit_berry_m2_cprior,fit_berry_m2_cprior2,fit_berry_m2_cprior3,
	QREJECTS,QREJECTS.CNT,
	qrejects0, qrejects0cp, qrejects0cp2, qrejects0cp3, qrejectsBerryM2, qrejectsBerryM2cp, qrejectsBerryM2cp2, qrejectsBerryM2cp3,
	pheno.sub, qqs_berry_m2, qqs_berry_m2cp, qqs_berry_m2cp2, qqs_berry_m2cp3, qqs0, qqs0cp, qqs0cp2, qqs0cp3, 
	cluster_info_tab, betas_sub, subAN, subBE, subMixBE, tmpBE, delta_mat, delta_cor_mat, delta_list
)
gc()


} ##########################################################################################################

