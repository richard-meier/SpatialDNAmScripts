
TM_CP2_3VARS = function(){
	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- mu[h,cli,j] + VX1[i]*gam1[h,cli,j] + VX2[i]*gam2[h,cli,j] + VX3[i]*gam3[h,cli,j]
					Z[i,h,cli,j] ~ dnorm(mmz[i,h,cli,j], psi[h,cli,j])
				}
				X[i,cli,j] ~ dnorm(mmm[i,cli,j],tau)
				mmm[i,cli,j] <- sum(WZ[i,,cli,j])
			}
		}
	}
	
	# =================== #
	# PRIOR SPECIFICATION #
	# =================== #
	
	tau ~ dgamma(0.844,0.001)
	for( cli in 1 : MCL ) {
		for( h in 1 : KCF ) {
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(0,0.01)
				gam1[h,cli,j] ~ dnorm(0,0.01)
				gam2[h,cli,j] ~ dnorm(0,0.01)
				gam3[h,cli,j] ~ dnorm(0,0.01)
				psi[h,cli,j] ~ dgamma(0.844,0.001)
			}
		}
	}
}

SCM2_CP2_3VARS = function(){
	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- DEGCL[cli]*omu[h,cli] + (1-DEGCL[cli])*mu[h,cli,j]
				wrapper.gam1[h,cli,j] <- DEGCL[cli]*ogam1[h,cli] + (1-DEGCL[cli])*gam1[h,cli,j]
				wrapper.gam2[h,cli,j] <- DEGCL[cli]*ogam2[h,cli] + (1-DEGCL[cli])*gam2[h,cli,j]
				wrapper.gam3[h,cli,j] <- DEGCL[cli]*ogam3[h,cli] + (1-DEGCL[cli])*gam3[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j] + VX1[i]*wrapper.gam1[h,cli,j] + VX2[i]*wrapper.gam2[h,cli,j] + VX3[i]*wrapper.gam3[h,cli,j]
					Z[i,h,cli,j] ~ dnorm(mmz[i,h,cli,j], psi[h,cli,j])
				}
				X[i,cli,j] ~ dnorm(mmm[i,cli,j],tau)
				mmm[i,cli,j] <- sum(WZ[i,,cli,j])
			}
		}
	}
	
	# =================== #
	# PRIOR SPECIFICATION #
	# =================== #
	
	tau ~ dgamma(0.844,0.001)
	otau ~ dgamma(0.1,0.1)
	ootau ~ dgamma(0.1,0.1)
	olambda1 ~ dgamma(0.1,0.1)
	olambda2 ~ dgamma(0.1,0.1)
	olambda3 ~ dgamma(0.1,0.1)
	oolambda1 ~ dgamma(0.1,0.1)
	oolambda2 ~ dgamma(0.1,0.1)
	oolambda3 ~ dgamma(0.1,0.1)
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(0,0.01)
		oogam1[h] ~ dnorm(0,0.01)
		oogam2[h] ~ dnorm(0,0.01)
		oogam3[h] ~ dnorm(0,0.01)
		for( cli in 1 : MCL ) {
			omu[h,cli] ~ dnorm(oomu[h],ootau)
			ogam1[h,cli] ~ dnorm(oogam1[h],oolambda1)
			ogam2[h,cli] ~ dnorm(oogam2[h],oolambda2)
			ogam3[h,cli] ~ dnorm(oogam3[h],oolambda3)
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(omu[h,cli],otau)
				gam1[h,cli,j] ~ dnorm(ogam1[h,cli],olambda1)
				gam2[h,cli,j] ~ dnorm(ogam2[h,cli],olambda2)
				gam3[h,cli,j] ~ dnorm(ogam3[h,cli],olambda3)
				psi[h,cli,j] ~ dgamma(0.844,0.001)
			}
		}
	}
}

generate.initial.values = function(clusters, chains, KCF=6, N=10){
	out = list()
	for(rrr in 1:chains){
		max_clst_length = max(unlist(lapply(clusters,FUN=length)))
		
		ooomu = runif(n=1, min=-0.5,max=0.5)
		otau = exp(runif(1,min=-1,max=1))
		ootau = exp(runif(1,min=-1,max=1))
		oootau = exp(runif(1,min=-1,max=1))
		
		ooogam1 = runif(n=1, min=-0.5,max=0.5)
		ooogam2 = runif(n=1, min=-0.5,max=0.5)
		ooogam3 = runif(n=1, min=-0.5,max=0.5)
		olambda1 = exp(runif(1,min=-1,max=1))
		olambda2 = exp(runif(1,min=-1,max=1))
		olambda3 = exp(runif(1,min=-1,max=1))
		oolambda1 = exp(runif(1,min=-1,max=1))
		oolambda2 = exp(runif(1,min=-1,max=1))
		oolambda3 = exp(runif(1,min=-1,max=1))
		ooolambda1 = exp(runif(1,min=-1,max=1))
		ooolambda2 = exp(runif(1,min=-1,max=1))
		ooolambda3 = exp(runif(1,min=-1,max=1))
		
		tau = exp(runif(1,min=1.5,max=3))
		
		mu = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		gam1 = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		gam2 = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		gam3 = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		psi = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		
		oomu = c()
		omu = array(NA,dim=c(KCF,length(clusters)))
		
		oogam1 = c()
		oogam2 = c()
		oogam3 = c()
		ogam1 = array(NA,dim=c(KCF,length(clusters)))
		ogam2 = array(NA,dim=c(KCF,length(clusters)))
		ogam3 = array(NA,dim=c(KCF,length(clusters)))
		
		Z = array(NA,dim=c(N,KCF,length(clusters),max_clst_length))
		
		for( h in 1 : KCF ) {
			oomu[h] = runif(n=1, min=-0.5,max=0.5)
			oogam1[h] = runif(n=1, min=-0.5,max=0.5)
			oogam2[h] = runif(n=1, min=-0.5,max=0.5)
			oogam3[h] = runif(n=1, min=-0.5,max=0.5)
		}
		for(ccc in 1:length(clusters)){
			for( h in 1 : KCF ) {
				omu[h,ccc] = runif(1,min=-0.1,0.1)
				ogam1[h,ccc] = runif(1,min=-0.1,0.1)
				ogam2[h,ccc] = runif(1,min=-0.1,0.1)
				ogam3[h,ccc] = runif(1,min=-0.1,0.1)
				for( j in 1 : length(clusters[[ccc]]) ) {
					mu[h,ccc,j] = runif(1,min=-0.1,0.1)
					gam1[h,ccc,j] = runif(1,min=-0.1,0.1)
					gam2[h,ccc,j] = runif(1,min=-0.1,0.1)
					gam3[h,ccc,j] = runif(1,min=-0.1,0.1)
					psi[h,ccc,j] = exp(runif(1,min=1.5,max=3))
					for(i in 1:N){
						Z[i,h,ccc,j] = runif(1,min=-0.1,0.1)
					}
				}
			}
		}
		clist = list(
			ooomu=ooomu, oomu=oomu, omu=omu, 
			ooogam1=ooogam1, ooogam2=ooogam2, ooogam3=ooogam3, 
			oogam1=oogam1, oogam2=oogam2, oogam3=oogam3, 
			ogam1=ogam1, ogam2=ogam2, ogam3=ogam3, 
			otau=otau, ootau=ootau, oootau=oootau, 
			olambda1=olambda1, olambda2=olambda2, olambda3=olambda3, 
			oolambda1=oolambda1, oolambda2=oolambda2, oolambda3=oolambda3, 
			ooolambda1=ooolambda1, ooolambda2=ooolambda2, ooolambda3=ooolambda3, 
			mu=mu, 
			gam1=gam1, gam2=gam2, gam3=gam3,
			psi=psi, tau=tau,
			Z=Z
		)
		out[[rrr]] = clist
	}
	return(out)
}

