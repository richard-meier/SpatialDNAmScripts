
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

calculate_super_cluster_span = function(scl){
	return( scl[[length(scl)]][length(scl[[length(scl)]])]-scl[[1]][1] )
}

draw_random_mixture_candidates = function(pheno){
	ctypes = c("Neutrophil", "NK", "Bcell", "CD4T", "CD8T","Monocyte")
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

prfun = function(x,a,b,c,minval,maxval){
	min_slope = b + 2*c*minval
	max_slope = b + 2*c*maxval
	min_y = a + b*minval + c*minval^2
	max_y = a + b*maxval + c*maxval^2
	min_extrap = min_y + (x-minval)*min_slope
	max_extrap = max_y + (x-maxval)*max_slope
	interp = a + b*x + c*x^2
	
	out_of_left_bound = x<minval
	out_of_right_bound = x>maxval
	within_bounds = !out_of_left_bound & !out_of_right_bound
	y = (out_of_left_bound)*min_extrap + (within_bounds)*interp + (out_of_right_bound)*max_extrap
	return(y)
}
