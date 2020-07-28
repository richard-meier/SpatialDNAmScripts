
generate.initial.values.d_old = function(clusters, chains, KCF=6, N=10){
	out = list()
	for(rrr in 1:chains){
		oomu = 0
		tomu = 1
		oogam = 0
		togam = 1
		
		tau = c()
		tmu = c()
		tgam = c()
		
		max_clst_length = max(unlist(lapply(clusters,FUN=length)))
		omu = array(NA,dim=c(length(clusters),KCF))
		ogam = array(NA,dim=c(length(clusters),KCF))
		
		mu = array(NA,dim=c(length(clusters),KCF,max_clst_length))
		gam = array(NA,dim=c(length(clusters),KCF,max_clst_length))
		psi = array(NA,dim=c(length(clusters),KCF,max_clst_length))
		
		XPRED = array(NA,dim=c(length(clusters),N,max_clst_length))
		
		Z = array(NA,dim=c(length(clusters),N,KCF,max_clst_length))
		
		for(ccc in 1:length(clusters)){
			tau[ccc] = exp(runif(1,min=-1,max=1))
			tmu[ccc] = exp(runif(1,min=-1,max=1))
			tgam[ccc] = exp(runif(1,min=-1,max=1))
			for(i in 1:N){
				for( j in 1 : length(clusters[[ccc]]) ) {
					XPRED[ccc,i,j] = runif(1,min=0.001,0.999)
				}
			}
			for( h in 1 : KCF ) {
				omu[ccc,h] = runif(1,min=-0.1,0.1)
				ogam[ccc,h] = runif(1,min=-0.1,0.1)
				for( j in 1 : length(clusters[[ccc]]) ) {
					mu[ccc,h,j] = runif(1,min=-0.1,0.1)
					gam[ccc,h,j] = runif(1,min=-0.1,0.1)
					psi[ccc,h,j] = exp(runif(1,min=-1,max=1))
					for(i in 1:N){
						Z[ccc,i,h,j] = runif(1,min=-0.1,0.1)
					}
				}
			}
		}
		clist = list(
			oomu=oomu, tomu=tomu, oogam=oogam, togam=togam,
			omu=omu, tmu=tmu, ogam=ogam, tgam=tgam,
			mu=mu, gam=gam, psi=psi,
			Z=Z, XPRED=XPRED, tau=tau
		)
		out[[rrr]] = clist
	}
	return(out)
}

generate.initial.values.berry.m1.old = function(clusters, chains, KCF=6, N=10){
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
		
		tau = exp(runif(1,min=-1,max=1))
		
		mu = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		gam = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		psi = array(NA,dim=c(length(clusters),KCF,max_clst_length))
		
		oomu = c()
		omu = array(NA,dim=c(KCF,length(clusters)))
		
		oogam = c()
		ogam = array(NA,dim=c(KCF,length(clusters)))
		
		XPRED = array(NA,dim=c(length(clusters),N,max_clst_length))
		QPRED = array(NA,dim=c(length(clusters),N,max_clst_length))
		
		Z = array(NA,dim=c(length(clusters),N,KCF,max_clst_length))
		QZ = array(NA,dim=c(length(clusters),N,KCF,max_clst_length))
		
		for( h in 1 : KCF ) {
			oomu[h] = runif(n=1, min=-0.5,max=0.5)
			oogam[h] = runif(n=1, min=-0.5,max=0.5)
		}
		for(ccc in 1:length(clusters)){
			for(i in 1:N){
				for( j in 1 : length(clusters[[ccc]]) ) {
					XPRED[ccc,i,j] = runif(1,min=0.001,0.999)
					QPRED[ccc,i,j] = runif(1,min=0.001,0.999)
				}
			}
			for( h in 1 : KCF ) {
				omu[h,ccc] = runif(1,min=-0.1,0.1)
				ogam[h,ccc] = runif(1,min=-0.1,0.1)
				for( j in 1 : length(clusters[[ccc]]) ) {
					mu[h,ccc,j] = runif(1,min=-0.1,0.1)
					gam[h,ccc,j] = runif(1,min=-0.1,0.1)
					psi[ccc,h,j] = exp(runif(1,min=-1,max=1))
					for(i in 1:N){
						Z[ccc,i,h,j] = runif(1,min=-0.1,0.1)
						QZ[ccc,i,h,j] = runif(1,min=-0.1,0.1)
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

generate.initial.values.berry.m1 = function(clusters, chains, KCF=6, N=10){
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

generate.initial.values.d = function(clusters, chains, KCF=6, N=10){
	out = list()
	for(rrr in 1:chains){
		oomu = 0
		tomu = 1
		oogam = 0
		togam = 1
		
		tau = exp(runif(1,min=1.5,max=3))
		tmu = c()
		tgam = c()
		
		max_clst_length = max(unlist(lapply(clusters,FUN=length)))
		omu = array(NA,dim=c(KCF,length(clusters)))
		ogam = array(NA,dim=c(KCF,length(clusters)))
		
		mu = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		gam = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		psi = array(NA,dim=c(KCF,length(clusters),max_clst_length))
		
		XPRED = array(NA,dim=c(N,length(clusters),max_clst_length))
		QPRED = array(NA,dim=c(N,length(clusters),max_clst_length))
		
		Z = array(NA,dim=c(N,KCF,length(clusters),max_clst_length))
		QZ = array(NA,dim=c(N,KCF,length(clusters),max_clst_length))
		
		for(ccc in 1:length(clusters)){
			tmu[ccc] = exp(runif(1,min=-1,max=1))
			tgam[ccc] = exp(runif(1,min=-1,max=1))
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
			oomu=oomu, tomu=tomu, oogam=oogam, togam=togam,
			omu=omu, tmu=tmu, ogam=ogam, tgam=tgam,
			mu=mu, gam=gam, psi=psi,
			Z=Z, XPRED=XPRED, tau=tau
		)
		out[[rrr]] = clist
	}
	return(out)
}
