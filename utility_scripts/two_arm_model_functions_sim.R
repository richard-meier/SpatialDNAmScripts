
TM_CP2 = function(){
	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- mu[h,cli,j] + GRP[i]*gam[h,cli,j]
					Z[i,h,cli,j] ~ dnorm(mmz[i,h,cli,j], psi[h,cli,j])
				}
				X[i,cli,j] ~ dnorm(mmm[i,cli,j],tau)
				mmm[i,cli,j] <- sum(WZ[i,,cli,j])
				
				XERR[i,cli,j] <- (X[i,cli,j] - XPRED[i,cli,j]) * (X[i,cli,j] - XPRED[i,cli,j])
				XPRED[i,cli,j] ~ dnorm(mmm[i,cli,j],tau)
			}
			XERRS1[i,cli] <- sum(XERR[i,cli,1:MCPG[cli]])
		}
		XERRS2[i] <- sum(XERRS1[i,])
	}
	SET <- sum(XERRS2[])
	MCPGSUM <- sum(MCPG[])
	RMSE <- sqrt(SET / (N*KCF*MCPGSUM))
	
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					qWZ[i,h,cli,j] <- QCFR[i,h] * QZ[i,h,cli,j]
					qmmz[i,h,cli,j] <- mu[h,cli,j] + GRP[i]*gam[h,cli,j]
					QZ[i,h,cli,j] ~ dnorm(qmmz[i,h,cli,j], psi[h,cli,j])
				}
				qmmm[i,cli,j] <- sum(qWZ[i,,cli,j])
				
				QERR[i,cli,j] <- (Q[i,cli,j] - QPRED[i,cli,j]) * (Q[i,cli,j] - QPRED[i,cli,j])
				QPRED[i,cli,j] ~ dnorm(qmmm[i,cli,j],tau)
			}
			QERRS1[i,cli] <- sum(QERR[i,cli,1:MCPG[cli]])
		}
		QERRS2[i] <- sum(QERRS1[i,])
	}
	
	QSET <- sum(QERRS2[])
	PREDRMSE <- sqrt(QSET / (N*KCF*MCPGSUM))
	
	# =================== #
	# PRIOR SPECIFICATION #
	# =================== #
	
	tau ~ dgamma(0.844,0.001)
	for( cli in 1 : MCL ) {
		for( h in 1 : KCF ) {
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(0,0.01)
				gam[h,cli,j] ~ dnorm(0,0.01)
				psi[h,cli,j] ~ dgamma(0.844,0.001)
			}
		}
	}
}

SCM2_CP2 = function(){
	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- DEGCL[cli]*omu[h,cli] + (1-DEGCL[cli])*mu[h,cli,j]
				wrapper.gam[h,cli,j] <- DEGCL[cli]*ogam[h,cli] + (1-DEGCL[cli])*gam[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j] + GRP[i]*wrapper.gam[h,cli,j]
					Z[i,h,cli,j] ~ dnorm(mmz[i,h,cli,j], psi[h,cli,j])
				}
				mmm[i,cli,j] <- sum(WZ[i,,cli,j])
				X[i,cli,j] ~ dnorm(mmm[i,cli,j],tau)
				XERR[i,cli,j] <- (X[i,cli,j] - XPRED[i,cli,j]) * (X[i,cli,j] - XPRED[i,cli,j])
				XPRED[i,cli,j] ~ dnorm(mmm[i,cli,j],tau)
			}
			XERRS1[i,cli] <- sum(XERR[i,cli,1:MCPG[cli]])
		}
		XERRS2[i] <- sum(XERRS1[i,])
	}
	SET <- sum(XERRS2[])
	MCPGSUM <- sum(MCPG[])
	RMSE <- sqrt(SET / (N*KCF*MCPGSUM))
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					qWZ[i,h,cli,j] <- QCFR[i,h] * QZ[i,h,cli,j]
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j] + GRP[i]*wrapper.gam[h,cli,j]
					QZ[i,h,cli,j] ~ dnorm(qmmz[i,h,cli,j], psi[h,cli,j])
				}
				qmmm[i,cli,j] <- sum(qWZ[i,,cli,j])
				QERR[i,cli,j] <- (Q[i,cli,j] - QPRED[i,cli,j]) * (Q[i,cli,j] - QPRED[i,cli,j])
				QPRED[i,cli,j] ~ dnorm(qmmm[i,cli,j],tau)
			}
			QERRS1[i,cli] <- sum(QERR[i,cli,1:MCPG[cli]])
		}
		QERRS2[i] <- sum(QERRS1[i,])
	}
	QSET <- sum(QERRS2[])
	PREDRMSE <- sqrt(QSET / (N*KCF*MCPGSUM))
	
	# =================== #
	# PRIOR SPECIFICATION #
	# =================== #
	
	tau ~ dgamma(0.844,0.001)
	otau ~ dgamma(0.1,0.1)
	ootau ~ dgamma(0.1,0.1)
	olambda ~ dgamma(0.1,0.1)
	oolambda ~ dgamma(0.1,0.1)
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(0,0.01)
		oogam[h] ~ dnorm(0,0.01)
		for( cli in 1 : MCL ) {
			omu[h,cli] ~ dnorm(oomu[h],ootau)
			ogam[h,cli] ~ dnorm(oogam[h],oolambda)
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(omu[h,cli],otau)
				gam[h,cli,j] ~ dnorm(ogam[h,cli],olambda)
				psi[h,cli,j] ~ dgamma(0.844,0.001)
			}
		}
	}
}

generate.initial.values.sim = function(clusters, chains, KCF=6, N=10){
	out = list()
	for(rrr in 1:chains){
		max_clst_length = max(unlist(lapply(clusters,FUN=length)))
		
		ooomu = runif(n=1, min=-0.5,max=0.5)
		otau = exp(runif(1,min=-1,max=1))
		ootau = exp(runif(1,min=-1,max=1))
		oootau = exp(runif(1,min=-1,max=1))
		ooogam = runif(n=1, min=-0.5,max=0.5)
		olambda = exp(runif(1,min=-1,max=1))
		oolambda = exp(runif(1,min=-1,max=1))
		ooolambda = exp(runif(1,min=-1,max=1))
		tau = exp(runif(1,min=1.5,max=3))
		
		mu = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		gam = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		psi = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		oomu = c()
		omu = array(NA,dim=c(KCF,length(clusters)))
		oogam = c()
		ogam = array(NA,dim=c(KCF,length(clusters)))
		XPRED = array(NA,dim=c(N,length(clusters),max_clst_length))
		QPRED = array(NA,dim=c(N,length(clusters),max_clst_length))
		Z = array(NA,dim=c(N,KCF,length(clusters),max_clst_length))
		QZ = array(NA,dim=c(N,KCF,length(clusters),max_clst_length))
		
		for( h in 1 : KCF ) {
			oomu[h] = runif(n=1, min=-0.5,max=0.5)
			oogam[h] = runif(n=1, min=-0.5,max=0.5)
		}
		for(ccc in 1:length(clusters)){
			for(i in 1:N){
				for( j in 1 : length(clusters[[ccc]]) ) {
					XPRED[i,ccc,j] = runif(1,min=0.001,0.999)
					QPRED[i,ccc,j] = runif(1,min=0.001,0.999)
				}
			}
			for( h in 1 : KCF ) {
				omu[h,ccc] = runif(1,min=-0.1,0.1)
				ogam[h,ccc] = runif(1,min=-0.1,0.1)
				for( j in 1 : length(clusters[[ccc]]) ) {
					mu[h,ccc,j] = runif(1,min=-0.1,0.1)
					gam[h,ccc,j] = runif(1,min=-0.1,0.1)
					psi[h,ccc,j] = exp(runif(1,min=1.5,max=3))
					for(i in 1:N){
						Z[i,h,ccc,j] = runif(1,min=-0.1,0.1)
						QZ[i,h,ccc,j] = runif(1,min=-0.1,0.1)
					}
				}
			}
		}
		clist = list(
			ooomu=ooomu, oomu=oomu, omu=omu, 
			ooogam=ooogam, oogam=oogam, ogam=ogam, 
			otau=otau, ootau=ootau, oootau=oootau, 
			olambda=olambda, oolambda=oolambda, ooolambda=ooolambda, 
			mu=mu, gam=gam, psi=psi, tau=tau,
			Z=Z, XPRED=XPRED, QZ=QZ, QPRED=QPRED
		)
		out[[rrr]] = clist
	}
	return(out)
}
