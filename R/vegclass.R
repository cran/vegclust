vegclass<-function(y, x) {
	if(!inherits(y,"vegclust")) stop("y must be a vegclust object")
	

	if(y$mode=="raw"){
		if(!is.null(y$fixedCenters)) centers = as.data.frame(rbind(y$mobileCenters,y$fixedCenters))
		else centers = as.data.frame(y$mobileCenters)
		if(length(names(x))!=length(names(centers)) || sum(names(x)==names(centers))<ncol(x)) {
				c = conformveg(x,centers)
				x = as.matrix(c$x)
				centers = as.matrix(c$y)
		} else {
			x = as.matrix(x)
			centers = as.matrix(centers)
		}
		k = nrow(centers)
	} else {
		memb = as.matrix(y$memb)
		m = y$m
		d2cl = as.matrix(y$dist2clusters)
		k = ncol(d2cl)
	}

	n = nrow(x)
		

   if(y$method=="NC") {
   	u = matrix(0,nrow=n,ncol=(k+1))
  	} else {
   	u = matrix(0,nrow=n,ncol=k)
  	}
	

	dist2cent = matrix(0,nrow=nrow(x),ncol=k)

	#1. compute distance to centroids (fixed and mobile)
	if(y$mode=="raw") {
	  	for(i in 1:k) {
		   dist2cent[,i] = sqrt(rowSums(sweep(x,2,centers[i,],"-")^2))
		}
	} else { #distance mode
		x = as.matrix(x)
		for(i in 1:k) {
			vargeom = sum((memb[,i]^m) %*% (d2cl[,i]^2))/sum(memb[,i]^m)
			for(j in 1:n) {
	  			dist2cent[j,i] = sqrt((sum((memb[,i]^m)*(x[j,]^2))/sum(memb[,i]^m))-vargeom)
			}
		}
	}
	
	  #2. compute membership to centroids
   	if (y$method=="KM") {
   		minv<-apply(dist2cent,1,min)
			for(k in 1:ncol(dist2cent)) u[,k] = as.numeric(dist2cent[,k]==minv)
			u[dist2cent==0]=1
		}
		else if(y$method=="NC") {
   		d2cm2<-cbind(dist2cent,y$dnoise)^2
   		for(k in 1:ncol(d2cm2)) {
   			a<-sweep(d2cm2,1,d2cm2[,k],"/")
   			u[,k] = 1/rowSums(a^(-1/(y$m-1)))
   		}
			u[d2cm2==0]=1
		} else if (y$method=="FCM") {
   		d2cm2<-dist2cent^2
   		for(k in 1:ncol(dist2cent)) {
   			a<-sweep(d2cm2,1,d2cm2[,k],"/")
   			u[,k] = 1/rowSums(a^(-1/(y$m-1)))
   		}
			u[dist2cent ==0]=1
		} else if (y$method=="PCM") {
			for(k in 1:ncol(dist2cent)) u[,k] = 1/(1+((dist2cent[,k]^2)/y$eta[k])^(1/(y$m-1)))
			u[dist2cent==0]=1
		} 
				
   #Prepare output
   u = as.data.frame(u)   
   dist2cent = as.data.frame(dist2cent)   
   if(ncol(u)==ncol(y$memb)) names(u) = names(y$memb)
   else {
   	names(u)[1:ncol(y$memb)] = names(y$memb)
   	names(u)[ncol(y$memb)+1] = "N"
   }
   names(dist2cent) = names(y$dist2clusters)   
	rownames(u) = rownames(x)
	rownames(dist2cent) = rownames(x)
	
   res = list(method = y$method, m =y$m, dnoise = y$dnoise, eta= y$eta, memb=u,dist2clusters=dist2cent)
   class(res)<-"vegclass"
	return(res)
		
}