
find.clusters.old = function(x,dth){
	out = list()
	out_idx = 1
	ctmp = c(x[1])
	for(i in 1:length(x)){
		if(x[i] > dth){
			out[[out_idx]] = ctmp
			ctmp = c(x[i])
			out_idx = out_idx + 1
		} else{
			ctmp = c(ctmp,x[i])
		}
	}
	out[[out_idx]] = ctmp
	return(out)
}

find.clusters = function(pos,dth){
	out = list()
	out_idx = 1
	ctmp = c(pos[1])
	for(i in 2:length(pos)){
		cdiff = pos[i]-pos[i-1]
		if(cdiff > dth){
			out[[out_idx]] = ctmp
			ctmp = c(pos[i])
			out_idx = out_idx + 1
		} else{
			ctmp = c(ctmp,pos[i])
		}
	}
	out[[out_idx]] = ctmp
	return(out)
}

split.cluster = function(cl,max_split_prop){
	distances = cl[-1] - cl[-length(cl)]
	distances = sort(distances,decreasing=TRUE)
	current_splits = list()
	for(i in 1:length(distances)){
		target_dist = distances[i]
		dist = 0
		k = 1
		while(dist != target_dist){
			k = k+1
			dist = cl[k] - cl[k-1]
		}
		current_splits = list(cl[1:(k-1)],cl[k:length(cl)])
		split_lengths = unlist(lapply(current_splits,FUN=length))
		if( max(split_lengths)/length(cl) <= max_split_prop ){
			break;
		}
	}
	return(current_splits)
}

split.large.clusters = function(clusters,max_cpg_num,max_split_prop){
	altered = FALSE
	clengths = unlist(lapply(clusters,FUN=length))
	new_clusts = list()
	k = 1
	for(i in 1:length(clusters)){
		to_be_added = list()
		if (clengths[i] > max_cpg_num){
			altered=TRUE
			to_be_added = split.cluster( clusters[[i]],max_split_prop )
		} else{
			to_be_added = list( clusters[[i]] )
		}
		for(r in 1:length(to_be_added)){
			new_clusts[[k]] = to_be_added[[r]]
			k = k+1
		}
	}
	return(list(clusters=new_clusts, was.split.performed=altered))
}

split.super.cluster = function(clusters,max_split_prop){
	distances = c()
	for(i in 2:length(clusters)){
		distances = c(distances, clusters[[i]][1] - clusters[[i-1]][length(clusters[[i-1]])] )
	}
	distances = sort(distances,decreasing=TRUE)	
	current_splits = list()
	for(i in 1:length(distances)){
		target_dist = distances[i]
		dist = 0
		k = 1
		while(dist != target_dist){
			k = k+1
			dist = clusters[[k]][1] - clusters[[k-1]][length(clusters[[k-1]])]
		}
		
		current_splits = list(clusters[1:(k-1)],clusters[k:length(clusters)])
		split_lengths = unlist(lapply(current_splits,FUN=length))
		if( max(split_lengths)/length(clusters) <= max_split_prop ){
			break;
		}
	}
	return(current_splits)
}

split.large.super.clusters = function(sclusts,max_cluster_num,max_split_prop){
	altered = FALSE
	clengths = unlist(lapply(sclusts,FUN=length))
	new_clusts = list()
	k = 1
	for(i in 1:length(sclusts)){
		to_be_added = list()
		if (clengths[i] > max_cluster_num){
			altered=TRUE
			to_be_added = split.super.cluster( sclusts[[i]],max_split_prop )
		} else{
			to_be_added = list( sclusts[[i]] )
		}
		for(r in 1:length(to_be_added)){
			new_clusts[[k]] = to_be_added[[r]]
			k = k+1
		}
	}
	return(list(clusters=new_clusts, was.split.performed=altered))
}

construct.cluster.structures = function(chr_positions,cl_threshold1,cl_threshold2,max_cn,max_scn,max_split_prop){
	clusters_naive = find.clusters(chr_positions,dth=cl_threshold1)
	clusters = clusters_naive
	split_performed = TRUE
	while(split_performed){
		csplit = split.large.clusters(clusters,max_cpg_num=max_cn,max_split_prop=max_split_prop)
		clusters = csplit$clusters
		split_performed = csplit$was.split.performed
	}
	
	supercl_v1 = find.super.clusters(clusters_naive,cl_threshold2)
	supercl_v2 = find.super.clusters(clusters,cl_threshold2)
	supercl_v3 = supercl_v2
	split_performed = TRUE
	while(split_performed){
		csplit = split.large.super.clusters(supercl_v3,max_cluster_num=max_scn,max_split_prop=max_split_prop)
		supercl_v3 = csplit$clusters
		split_performed = csplit$was.split.performed
	}
	
	return( list(sc.v1=supercl_v1,sc.v2=supercl_v2,sc.v3=supercl_v3) )
}

find.super.clusters = function(clusters,dth){
	out = list()
	out_idx = 1
	ctmp = list(clusters[[1]])
	for(i in 2:length(clusters)){
		cpos1 = tail(clusters[[i-1]],n=1)
		cpos2 = head(clusters[[i]],n=1)
		cdiff = cpos2-cpos1
		if(cdiff > dth){
			out[[out_idx]] = ctmp
			ctmp = list(clusters[[i]])
			out_idx = out_idx + 1
		} else{
			ctmp = c(ctmp,clusters[i])
		}
	}
	out[[out_idx]] = ctmp
	return(out)
}

find.cluster.length.old = function(x,dth){
	out = c()
	out_idx = 1
	clen = 1
	for(i in 1:length(x)){
		if(x[i] > dth){
			out[out_idx] = clen
			clen = 1
			out_idx = out_idx + 1
		} else{
			clen = clen + 1
		}
	}
	out[out_idx] = clen
	return(out)
}

calculate_super_cluster_span = function(scl){
	return( scl[[length(scl)]][length(scl[[length(scl)]])]-scl[[1]][1] )
}

draw_random_mixture_candidates = function(pheno){
	ctypes = unique(pheno$Cell.Type)
	indices = c()
	for(ctp in ctypes){
		indices = c(indices, sample(which(pheno$Cell.Type == ctp))[1])
	}
	return(indices)
}

sigmoid = function(x){
	return( 1/(1+exp(-x) ) )
}
logit = function(x){
	return( log(x) - log(1-x) )
}

llt.private = function(x, delta=0, noise = 0.1){
	return( sigmoid(logit(x) + rnorm(n=1,mean=delta,sd=noise)) )
}

logit.linear.transform = function(x, delta=0, noise = 0.1){
	if(length(x)>1){
		return( sapply(x,FUN=llt.private, delta=delta, noise=noise) )
	} else{
		return( llt.private(x=x, delta=delta, noise=noise) )
	}
}

check.larger.0 = function(x){
	return( sum(x>0) == length(x) )
}

autocor = function(lag,chain){
	if(lag<=0){
		return(1)
	} else{
		ch1 = chain[-(1:lag)]
		ch2 = chain[-((length(chain)-lag+1):length(chain))]
		return(cor(ch1,ch2))
	}
}

cross.rate = function(chain,normalized=TRUE){
	fff = 0.25
	qts = quantile(chain, probs=c(fff,1-fff))
	lower = FALSE
	upper = FALSE
	cross = 0
	for(j in 1:length(chain)){
		if(chain[j] < qts[1]){
			if(!lower){
				cross = cross+1
			}
			lower = TRUE
		} else if(chain[j] > qts[2]){
			if(!upper){
				cross = cross+1
			}
			upper = TRUE
		} else{
			lower = FALSE
			upper = FALSE
		}
	}
	if(normalized){
		return( 3*cross / length(chain) )
	} else{
		return( cross / length(chain) )
	}
}

create.combinations = function(...){
	x <- list(...)
	return(ccbn(x,1))
}
ccbn = function(xlist,n){
	if(n >= length(xlist)){
		out = data.frame(xlist[[length(xlist)]])
		colnames(out) = names(xlist)[length(xlist)]
		return( out )
	} else{
		base = ccbn(xlist,n+1)
		out = data.frame(xlist[[n]][1], base)
		colnames(out) = c( names(xlist)[n], colnames(base) )
		if( length(xlist[[n]])>1 ){
			for(i in 2:length(xlist[[n]])){
				tmp = data.frame(xlist[[n]][i], base)
				colnames(tmp) = c( names(xlist)[n], colnames(base) )
				out = rbind(out,tmp)
			}
		}
		return(out)
	}
}

diag.line=function(col,lty=1,identity=TRUE){
	if(identity){
		lines(c(-10^10,10^10),c(-10^10,10^10),col=col,lty=lty)
	} else{
		lines(c(-10^10,10^10),c(10^10,-10^10),col=col,lty=lty)
	}
}
