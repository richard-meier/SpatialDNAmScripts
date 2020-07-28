
marginal_berry_like_model_m3_custom_prior2 = function(){

	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- mu[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
	
	ootau ~ dgamma(0.1,0.1)
	
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(0,0.01)
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(oomu[h],ootau)
				psi[h,cli,j] ~ dgamma(0.844,0.001)
			}
		}
	}
	
	tau ~ dgamma(0.844,0.001)
}



marginal_berry_like_model_m2_custom_prior2 = function(){

	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- DEGCL[cli]*omu[h,cli] + (1-DEGCL[cli])*mu[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
	
	otau ~ dgamma(0.1,0.1)
	ootau ~ dgamma(0.1,0.1)
	
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(0,0.01)
		for( cli in 1 : MCL ) {
			omu[h,cli] ~ dnorm(oomu[h],ootau)
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(omu[h,cli],otau)
				psi[h,cli,j] ~ dgamma(0.844,0.001)
			}
		}
	}
	
	tau ~ dgamma(0.844,0.001)
}



marginal_berry_like_model_m1_custom_prior2 = function(){

	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- DEGCL[cli]*omu[h,cli] + (1-DEGCL[cli])*mu[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
	
	ooomu ~ dnorm(0,0.01)
	otau ~ dgamma(0.1,0.1)
	ootau ~ dgamma(0.1,0.1)
	oootau ~ dgamma(0.1,0.1)
	
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(ooomu,oootau)
		for( cli in 1 : MCL ) {
			omu[h,cli] ~ dnorm(oomu[h],ootau)
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(omu[h,cli],otau)
				psi[h,cli,j] ~ dgamma(0.844,0.001)
			}
		}
	}
	
	tau ~ dgamma(0.844,0.001)
}



marginal_tca_like_model_flip_custom_prior2 = function(){
	
	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- mu[h,cli,j]
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
				psi[h,cli,j] ~ dgamma(0.844,0.001)
			}
		}
	}
}



marginal_berry_like_model_m3_custom_prior = function(){

	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- mu[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
	
	ootau ~ dgamma(0.1,0.1)
	
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(0,0.01)
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(oomu[h],ootau)
				psi[h,cli,j] ~ dgamma(1.36,0.01)
			}
		}
	}
	
	tau ~ dgamma(1.36,0.01)
}



marginal_berry_like_model_m2_custom_prior = function(){

	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- DEGCL[cli]*omu[h,cli] + (1-DEGCL[cli])*mu[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
	
	otau ~ dgamma(0.1,0.1)
	ootau ~ dgamma(0.1,0.1)
	
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(0,0.01)
		for( cli in 1 : MCL ) {
			omu[h,cli] ~ dnorm(oomu[h],ootau)
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(omu[h,cli],otau)
				psi[h,cli,j] ~ dgamma(1.36,0.01)
			}
		}
	}
	
	tau ~ dgamma(1.36,0.01)
}



marginal_berry_like_model_m1_custom_prior = function(){

	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- DEGCL[cli]*omu[h,cli] + (1-DEGCL[cli])*mu[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
	
	ooomu ~ dnorm(0,0.01)
	otau ~ dgamma(0.1,0.1)
	ootau ~ dgamma(0.1,0.1)
	oootau ~ dgamma(0.1,0.1)
	
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(ooomu,oootau)
		for( cli in 1 : MCL ) {
			omu[h,cli] ~ dnorm(oomu[h],ootau)
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(omu[h,cli],otau)
				psi[h,cli,j] ~ dgamma(1.36,0.01)
			}
		}
	}
	
	tau ~ dgamma(1.36,0.01)
}



marginal_tca_like_model_flip_custom_prior = function(){
	
	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- mu[h,cli,j]
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
	
	tau ~ dgamma(1.36,0.01)
	for( cli in 1 : MCL ) {
		for( h in 1 : KCF ) {
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(0,0.01)
				psi[h,cli,j] ~ dgamma(1.36,0.01)
			}
		}
	}
}



marginal_berry_like_model_m3 = function(){

	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- mu[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
	
	ootau ~ dgamma(0.1,0.1)
	
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(0,0.01)
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(oomu[h],ootau)
				psi[h,cli,j] ~ dgamma(0.01,0.01)
			}
		}
	}
	
	tau ~ dgamma(0.01,0.01)
}



marginal_berry_like_model_m2 = function(){

	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- DEGCL[cli]*omu[h,cli] + (1-DEGCL[cli])*mu[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
	
	otau ~ dgamma(0.1,0.1)
	ootau ~ dgamma(0.1,0.1)
	
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(0,0.01)
		for( cli in 1 : MCL ) {
			omu[h,cli] ~ dnorm(oomu[h],ootau)
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(omu[h,cli],otau)
				psi[h,cli,j] ~ dgamma(0.01,0.01)
			}
		}
	}
	
	tau ~ dgamma(0.01,0.01)
}



marginal_berry_like_model_m1 = function(){

	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	
	for( cli in 1 : MCL ) {
		for( j in 1 : MCPG[cli] ) {
			for( h in 1 : KCF ) {
				wrapper.mu[h,cli,j] <- DEGCL[cli]*omu[h,cli] + (1-DEGCL[cli])*mu[h,cli,j]
			}
		}
	}
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- wrapper.mu[h,cli,j]
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
	
	ooomu ~ dnorm(0,0.01)
	otau ~ dgamma(0.1,0.1)
	ootau ~ dgamma(0.1,0.1)
	oootau ~ dgamma(0.1,0.1)
	
	for( h in 1 : KCF ) {
		oomu[h] ~ dnorm(ooomu,oootau)
		for( cli in 1 : MCL ) {
			omu[h,cli] ~ dnorm(oomu[h],ootau)
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(omu[h,cli],otau)
				psi[h,cli,j] ~ dgamma(0.01,0.01)
			}
		}
	}
	
	tau ~ dgamma(0.01,0.01)
}



marginal_tca_like_model_flip = function(){
	
	# ======================== #
	# LIKELIHOOD SPECIFICATION #
	# ======================== #
	for( i in 1 : N ) {
		for( cli in 1 : MCL ) {
			for( j in 1 : MCPG[cli] ) {
				for( h in 1 : KCF ) {
					WZ[i,h,cli,j] <- CFR[i,h] * Z[i,h,cli,j]
					mmz[i,h,cli,j] <- mu[h,cli,j]
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
					qmmz[i,h,cli,j] <- mu[h,cli,j]
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
	
	tau ~ dgamma(0.01,0.01)
	for( cli in 1 : MCL ) {
		for( h in 1 : KCF ) {
			for( j in 1 : MCPG[cli] ) {
				mu[h,cli,j] ~ dnorm(0,0.01)
				psi[h,cli,j] ~ dgamma(0.01,0.01)
			}
		}
	}
}
