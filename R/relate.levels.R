relate.levels<-function(lower, upper, verbose=FALSE, ...) {
  minupperclasses =999999
  maxupperclasses =0
	minlowerclasses=999999
  	maxlowerclasses =0
	for(i in 1:length(upper)) {
		if(!is.null(upper[[i]])) {
			nclasses = length(upper[[i]]$size)
			maxupperclasses = max(maxupperclasses,nclasses)
			minupperclasses = min(minupperclasses,nclasses)
		}
	}
	for(i in 1:length(lower)) {
		if(!is.null(lower[[i]])) {
			nclasses = length(lower[[i]]$size)
			minlowerclasses = min(minlowerclasses,nclasses)
			maxlowerclasses = max(maxlowerclasses,nclasses)
		}
	}
	#Prepare output matrices
	nnoi = data.frame(matrix(NA,length(minupperclasses:maxupperclasses), length(minlowerclasses:maxlowerclasses)))
	row.names(nnoi) = minupperclasses:maxupperclasses
	names(nnoi) = minlowerclasses:maxlowerclasses
	miss = nnoi
	nint = nnoi
	missfix = nnoi
	missmob = nnoi
	
	#Loop over lower clustering
	for(j in 1:length(lower)) {
		if(!is.null(lower[[j]])) {
			nlowerclasses = length(lower[[j]]$size)
			if(verbose) cat(paste("Number of lower classes:", nlowerclasses,"\n"))
			if(!is.null(lower[[j]]$fixedCenters)) {
					centroids = data.frame(rbind(lower[[j]]$mobileCenters,lower[[j]]$fixedCenters))
			} else centroids = lower[[j]]$mobileCenters
			#Loop over upper clustering
			for(i in 1:length(upper)) {
				if(!is.null(upper[[i]])){
					if(verbose) cat(".")
					nupperclasses = length(upper[[i]]$size)
					#Assign lower centroids to upper classes and defuzzify
					memb <- defuzzify(vegclass(upper[[i]],centroids), ...)$memb
					rowind<-nupperclasses-minupperclasses+1
					colind<-nlowerclasses-minlowerclasses+1
					nint[rowind,colind]=sum(rowSums(memb)==0)
					miss[rowind,colind]=sum(colSums(memb[,1:nupperclasses])==0)
					if(upper[[i]]$method=="NC") nnoi[rowind,colind]=sum(memb[,ncol(memb)])
					if(!is.null(upper[[i]]$fixedCenters)) {
						nmob = nrow(upper[[i]]$mobileCenters)
						nfix = nrow(upper[[i]]$fixedCenters)
						if(nmob==1) missmob[rowind,colind]=as.numeric(sum(memb[,1])==0)
						else missmob[rowind,colind]=sum(colSums(memb[,1:nmob])==0)
						missfix[rowind,colind]=sum(colSums(memb[,(nmob+1):(nmob+nfix)])==0)
					}
				}
			}
		}
		if(verbose) cat("done.\n")
	}	

	return(list(noise=nnoi, inter=nint, miss=miss, missfix=missfix, missmob= missmob))
}