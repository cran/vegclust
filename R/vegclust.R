vegclust <-
function(x,mobileCenters, fixedCenters=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100, nstart=1, maxminJ=10, seeds = NULL, verbose=FALSE) {

#One run of vegclust   
vegclustone <-
function(x,mobileCenters, fixedCenters=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100) {
   METHODS <- c("KM", "FCM", "PCM","NC")
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
   
	x = as.matrix(x)
	n = nrow(x)
	
	#Sets the starting mobile centroids
	if(is.data.frame(mobileCenters) || is.matrix(mobileCenters)) {
		mobileCenters = as.matrix(mobileCenters)
		if(ncol(mobileCenters)!=ncol(x)) {
			stop("The number and identity of species for mobile centers must be the same as the number and identity of species for x")
		}		
	}
	else if(is.vector(mobileCenters) && length(mobileCenters)==1 && is.numeric(mobileCenters)) {
		if(mobileCenters==1) mobileCenters = t(as.matrix(x[sample(n,mobileCenters),]))
		else mobileCenters = as.matrix(x[sample(n,mobileCenters),])
	}
	else if(is.vector(mobileCenters) && is.numeric(mobileCenters)) {
		mobileCenters = as.matrix(x[mobileCenters,])
	} 
	if(!is.matrix(mobileCenters)) {
		stop("Provide a number, a vector of seeds, or coordinates for mobile centers")
	}	
	kMov = nrow(mobileCenters)
	rownames(mobileCenters)<-c(1:kMov)
	
	#Sets the fixed centroids
	if(!is.null(fixedCenters)) {
		if(is.data.frame(fixedCenters)) {
			fixedCenters = as.matrix(fixedCenters)
			kFix = nrow(fixedCenters)
		}
		else if(!is.matrix(fixedCenters)) {
			stop("Fixed centers must be specified as a matrix or a data frame")
		}	
		else {kFix = nrow(fixedCenters)}
		if(ncol(fixedCenters)!=ncol(x)) {
			stop("The number and identity of species for fixed centers must be the same as the number and identity of species for x")
		}
	} else {
		kFix = 0
	}

   if(method=="PCM" && length(eta)!=(kMov+kFix)) stop("Vector of reference distances (eta) must have a length equal to the number of clusters")
  	dist2cent = matrix(0,nrow=n,ncol=(kMov+kFix))

   if(method=="NC") {
   	u = matrix(0,nrow=n,ncol=(kMov+kFix+1))
  		uPrev = matrix(0,nrow=n,ncol=(kMov+kFix+1))
  	} else {
   	u = matrix(0,nrow=n,ncol=(kMov+kFix))
  		uPrev = matrix(0,nrow=n,ncol=(kMov+kFix))
  	}

	continue = TRUE
	iter = 1
   #iterates until no change in memberships
   while(continue) {
		#1. compute distance to centroids (fixed and mobile)
   	for(k in 1:kMov) {
   		dist2cent[,k] = sqrt(rowSums(sweep(x,2,mobileCenters[k,],"-")^2))
   	}
  		if(kFix==1) {
  			dist2cent[,1+kMov] = sqrt(rowSums(sweep(x,2, as.vector(fixedCenters),"-")^2))
  		}  			
  		else if(kFix>1) {
  			for(k in 1:kFix) {
  				dist2cent[,k+kMov] = sqrt(rowSums(sweep(x,2,fixedCenters[k,],"-")^2))
  			}
  		}  			
	  #2. compute membership to centroids
   	if (method=="KM") {
   		minv<-apply(dist2cent,1,min)
			for(k in 1:ncol(dist2cent)) u[,k] = as.numeric(dist2cent[,k]==minv)
			u[dist2cent==0]=1
		}
		else if(method=="NC") {
   		d2cm2<-cbind(dist2cent,dnoise)^2
   		for(k in 1:ncol(d2cm2)) {
   			a<-sweep(d2cm2,1,d2cm2[,k],"/")
   			u[,k] = 1/rowSums(a^(-1/(m-1)))
   		}
			u[d2cm2==0]=1
		} else if (method=="FCM") {
   		d2cm2<-dist2cent^2
   		for(k in 1:ncol(dist2cent)) {
   			a<-sweep(d2cm2,1,d2cm2[,k],"/")
   			u[,k] = 1/rowSums(a^(-1/(m-1)))
   		}
			u[dist2cent ==0]=1
		} else if (method=="PCM") {
			for(k in 1:ncol(dist2cent)) u[,k] = 1/(1+((dist2cent[,k]^2)/eta[k])^(1/(m-1)))
			u[dist2cent==0]=1
		} 	 	
   	#Check for stopping
   	if(iter>2) {
   		continue = (max(abs(u-uPrev))>alpha) && (iter<=iter.max) && (max(abs(u-uPrev2))>alpha)
   	}
   	
   	if(continue) {
	   	#3. update mobile centroids
	   	for(k in 1:kMov) mobileCenters[k,]=((u[,k]^m)%*%x)/sum(u[,k]^m)
	   	uPrev2 = uPrev
	   	uPrev = u	   	
	   	iter=iter+1
	   	if(verbose) cat(".")
   	}
   }
   if(verbose) cat(paste("\nIterations:", iter,"\n"))
  	if(method=="FCM" || method=="KM") functional = sum((dist2cent^2)*(u^m))
  	else if(method=="NC") functional = sum((dist2cent^2)*(u[,-(kMov+kFix+1)]^m))+sum(dnoise^2*u[,kMov+kFix+1]^m)
  	else if(method=="PCM") {
  		functional = 0
  		for(k in 1:(kMov+kFix)) functional = functional+sum((dist2cent[,k]^2)*(u[,k]^m))+sum(eta[k]*(1-u[,k])^m)
  	}
   
   #Prepare output
   mobileCenters = as.data.frame(mobileCenters)
   u = as.data.frame(u)   
   dist2cent = as.data.frame(dist2cent)   
	for(k in 1:kMov) {
		row.names(mobileCenters)[k] = paste("M",k,sep="")
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
	withinss = colSums((dist2cent^2)*(u[,1:(kMov+kFix)]^m))
   res = list(mode="raw", method=method, m = m, dnoise = dnoise,eta = eta, memb=u,mobileCenters=mobileCenters, fixedCenters=fixedCenters, dist2clusters=dist2cent, withinss = withinss, size=size, functional=functional, iter=iter)
   class(res)<-"vegclust"
	return(res)
}

	if(is.null(seeds)) seeds = 1:nrow(x)
	#print(seeds)
   #If mobileCenters is a number and nstart>1 perform different random starts
	if(is.vector(mobileCenters) && length(mobileCenters)==1 && is.numeric(mobileCenters)) {
	   bestRun = vegclustone(x,mobileCenters=x[sample(seeds,mobileCenters),], fixedCenters, method, m,dnoise, eta, alpha, iter.max)
		if(nstart>1) {
			minJ = 0
			i = 2
			while(i<=nstart) {
				run = vegclustone(x,mobileCenters=x[sample(seeds,mobileCenters),], fixedCenters, method, m,dnoise,eta, alpha, iter.max)
				if(run$functional<=bestRun$functional) {
					bestRun = run
					if(max(abs(run$memb-bestRun$memb))<alpha) {
						minJ = minJ+1
					} else {
						minJ = 1
					}
				}
				if(minJ==maxminJ) {
					i=nstart
					if(verbose) print("Maximum number of minimum J reached. Stopping.")
				}
				i=i+1
			}
		}		
		return(bestRun)
	} else { #Perform a single run
		return(vegclustone(x,mobileCenters, fixedCenters, method, m,dnoise, eta, alpha, iter.max))
	}
}

