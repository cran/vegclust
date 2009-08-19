vegclustdist <-
function(x,mobileMemb, fixedMemb=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100, nstart=1, seeds = NULL) {

#One run of vegclustdist   
vegclustonedist <-
function(d,mobileMemb, fixedMemb=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100) {
   METHODS <- c("KM", "FCM", "PCM", "NC")
   method <- match.arg(method, METHODS)
   if(method=="KM") {
   	m=1.0
   	dnoise=NULL
   	eta=NULL
   }
   else if(method=="FCM") {
   	dnoise=NULL
   	eta=NULL
   }
   else if(method=="NC") {
   	if(is.null(dnoise)) stop("Must provide a value for dnoise")
   	eta = NULL
   }
   else if(method=="PCM") {
   	if(is.null(eta)) stop("Must provide a vector of values for eta")
   	dnoise = NULL
   }
   
	d = as.matrix(d)
	n = nrow(d)
	
	#Sets the starting memberships for mobile clusters
	if(is.data.frame(mobileMemb) || is.matrix(mobileMemb)) {
		if(nrow(mobileMemb)!=ncol(d)) {
			stop("The number of rows in mobileMemb must be the same as the number rows and columns of d")
		}		
		u = as.matrix(mobileMemb)
	}
	else if(is.vector(mobileMemb) && length(mobileMemb)==1 && is.numeric(mobileMemb)) {
		s = sample(n, mobileMemb)
		u = matrix(0,n,length(s))
		for(k in 1:length(s)) {
			u[s[k],k]=1
		}
	} else if(is.vector(mobileMemb) && is.numeric(mobileMemb)) {
		u = matrix(0,n,length(mobileMemb))
		for(k in 1:length(mobileMemb)) {
			u[mobileMemb[k],k]=1
		}		
	}
	else {
		stop("Provide a number, a vector of seeds, or membership matrix for mobile clusters")
	}	
	kMov = ncol(u)
	
	#Sets the fixed cluster memberships
	if(!is.null(fixedMemb)) {
		if(is.data.frame(fixedMemb)) {
			kFix = ncol(fixedMemb)
			u = cbind(u,fixedMemb)
		}
		else if(!is.matrix(fixedMemb)) {
			stop("Fixed clusters must be specified as a matrix or a data frame")
		}	
		else {
			kFix = ncol(fixedMemb)
			u = cbind(u,fixedMemb)
		}
	} else {
		kFix = 0
	}
	#Check parameters
   if(method=="PCM" && length(eta)!=(kMov+kFix)) stop("Vector of reference distances (eta) must have a length equal to the number of clusters")

   #add extra column for NC
   if(method=="NC") {
   	u = cbind(u, vector("numeric",length=n))
  	}
  	uPrev = matrix(0,nrow=n,ncol=ncol(u))
	
   #initialize squared distances to cluster and sets the fixed squared distances
   sqdist2cent = matrix(0,nrow=n,ncol=(kMov+kFix))
   if(kFix>0) {
   	for(k in (kMov+1):(kMov+kFix)) {
			vargeom = sum((u[,k]^m) %*% (d^2) %*% (u[,k]^m))/(2*sum(u[,k]^m)^2)
	  		for(i in 1:n) {
	  			sqdist2cent[i,k] = (sum((u[,k]^m)*(d[i,]^2))/sum(u[,k]^m))-vargeom
	  			if(sqdist2cent[i,k]<0) sqdist2cent[i,k]=0
	  		}
   	}   	
   }

	continue = TRUE
	iter = 1
   #iterates until no change in memberships
   while(continue) {
   	vargeom = vector("numeric", kMov)
   	for(k in 1:(kMov)) {
			vargeom[k] = sum((u[,k]^m) %*% (d^2) %*% (u[,k]^m))/(2*sum(u[,k]^m)^2)
   	}
   	#for all objects
  		for(i in 1:n) {
		   #1. compute squared distances to mobile clusters
  			for(k in 1:kMov) {
	  			sqdist2cent[i,k] = (sum((u[,k]^m)*(d[i,]^2))/sum(u[,k]^m))-vargeom[k]
	  			if(sqdist2cent[i,k]<0) sqdist2cent[i,k]=0
  			}
	   	#2. compute membership to centroids for mobile and fixed clusters
			for(k in 1:(kMov+kFix)) {
				if(method=="NC") {
					u[i,k] = 1/(sum(((sqdist2cent[i,k])/sqdist2cent[i,])^(1/(m-1)))+ ((sqdist2cent[i,k])/dnoise^2)^(1/(m-1)))
					if(sqdist2cent[i,k]==0) u[i,k]=1
				} else if (method=="FCM") {
					u[i,k] = 1/(sum(((sqdist2cent[i,k])/sqdist2cent[i,])^(1/(m-1))))
					if(sqdist2cent[i,k]==0) u[i,k]=1
				} else if (method=="PCM") {
					u[i,k] = 1/(1+((sqdist2cent[i,k])/eta[k])^(1/(m-1)))
					if(sqdist2cent[i,k]==0) u[i,k]=1
				} else if (method=="KM") {
					u[i,k] = ifelse(sqdist2cent[i,k]==min(sqdist2cent[i,]),1,0)
				}
			}
			if(method=="NC") u[i,(kMov+kFix+1)] = 1-sum(u[i,1:(kMov+kFix)])
  	 	}
	  	 	
   	#Check for stopping
   	if(iter>2) {
   		continue = (max(abs(u-uPrev))>alpha) && (iter<=iter.max) && (max(abs(u-uPrev2))>alpha)
   	}
   	
   	if(continue) {
	   	uPrev2 = uPrev
	   	uPrev = u	   	
	   	iter=iter+1
   	}
   }
  	if(method=="FCM" || method=="KM") functional = sum((sqdist2cent)*(u^m))
  	else if(method=="NC") functional = sum((sqdist2cent)*(u[,-(kMov+kFix+1)]^m))+sum(dnoise^2*u[,kMov+kFix+1]^m)
  	else if(method=="PCM") {
  		functional = 0
  		for(k in 1:(kMov+kFix)) functional = functional+sum((sqdist2cent[,k])*(u[,k]^m))+sum(eta[k]*(1-u[,k])^m)
  	}
   
   #Prepare output
   u = as.data.frame(u)   
   dist2cent = as.data.frame(sqrt(sqdist2cent))   
	for(k in 1:kMov) {
		names(u)[k] = paste("M",k,sep="")
		names(dist2cent)[k] = paste("M",k,sep="")
	}
	if(kFix>1) {
		for(k in (kMov+1):(kMov+kFix)) {
			names(u)[k] = paste("F",k,sep="")
			names(dist2cent)[k] = paste("F",k,sep="")
		}
	}
	if(method=="NC") names(u)[kMov+kFix+1] = "N"
	rownames(u) = rownames(x)
	rownames(dist2cent) = rownames(x)
	size = colSums(u[,1:(kMov+kFix)])
	withinss = colSums((sqdist2cent)*(u[,1:(kMov+kFix)]^m))
   res = list(mode="dist", method=method, m = m, dnoise = dnoise,eta = eta, memb=u,mobileCenters=NULL, fixedCenters=NULL, dist2clusters=dist2cent, withinss = withinss, size=size, functional=functional)
   class(res)<-"vegclust"
	return(res)
}

	if(is.null(seeds)) seeds = 1:nrow(as.matrix(x))
   #If mobileCenters is a number and nstart>1 perform different random starts
	if(is.vector(mobileMemb) && length(mobileMemb)==1 && is.numeric(mobileMemb)) {
	   bestRun = vegclustonedist(x, mobileMemb=sample(seeds, mobileMemb), fixedMemb, method, m,dnoise, eta, alpha, iter.max)
		if(nstart>1) {
			for(i in 2:nstart) {
				run = vegclustonedist(x,mobileMemb=sample(seeds,mobileMemb), fixedMemb, method, m,dnoise,eta, alpha, iter.max)
				if(run$functional<bestRun$functional) {
					bestRun = run
				}
			}
		}		
		return(bestRun)
	} else { #Perform a single run
		return(vegclustonedist(x,mobileMemb, fixedMemb, method, m,dnoise, eta, alpha, iter.max))
	}
}

