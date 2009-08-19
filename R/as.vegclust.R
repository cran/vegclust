as.vegclust <-
function(x,cluster) {
   cln =levels(as.factor(cluster))
   if(inherits(x,"dist")) {
   	mode="distance"
   	x = as.matrix(x)
   } else {
   	mode="raw"
   }
   k = length(cln)
   n = nrow(x)
   
   dist2cent = matrix(0,nrow=n,ncol=k) 
   if(mode=="raw") centers = matrix(0,nrow=k,ncol=ncol(x)) 
   else centers = NULL
   
   dist2onecluster<-function(x,object) {
		#x is an (euclidean) distance matrix
		x = as.matrix(x)
		n = nrow(x)
		if (length(object)!=n) stop("Length of object vector must be equal to the number of sites in x")
		vargeom = (rep(1,n) %*% (x^2) %*% rep(1,n))/(2*(n^2))
		return(sqrt((sum(object^2)/n)-vargeom))
	}
   dist2clusters<-function(x,cluster,object) {
	   n = nrow(as.matrix(x))
		if (length(cluster)!=n) 
            stop("Length of cluster vector must be equal to the number of sites in x")
		cluster = as.factor(cluster)
		k = length(levels(cluster))
		d = vector("numeric",k)
		for(i in 1:k) {
			sel = (cluster==levels(cluster)[i])
		  sel[is.na(sel)]=FALSE
			d[i] = dist2onecluster(as.dist(as.matrix(x)[sel,sel]),object[sel])
		}	
		names(d) = levels(cluster)
		return(d)
   }	
	
   disteu <-function(x,y) {
   	return(sqrt(sum((x-y)^2)))
   }

   #Memberships
   u = data.frame(as.memb(cluster))
	 rownames(u) = rownames(x)
	 names(u) = levels(as.factor(cluster))
   
   #Distance to clusters
   if(mode=="distance") {
   	for(j in 1:n) {
   		dist2cent[j,] = dist2clusters(x,cluster,x[j,])
   	}
   } else{
   	for(i in 1:k) {
   		sel = (cluster==cln[i])
   		sel[is.na(sel)]=FALSE
   		centers[i,] = apply(x[sel,],2,"mean")
		  dist2cent[,i] = apply(x,1,"disteu",centers[i,])
   	}   
   }   
   
   #Centers
   if(mode=="raw"){
   	centers = as.data.frame(centers)   
   	names(centers) = names(x)
   	row.names(centers) = levels(as.factor(cluster))
   } 
   
   dist2cent = as.data.frame(dist2cent)  
	names(dist2cent) = levels(as.factor(cluster))
	rownames(dist2cent) = rownames(x)
	
	size = colSums(u)
	withinss = colSums((dist2cent^2)*u)
	functional = sum(withinss)
	
   res = list(mode = mode, method="KM", m = 1.0, dnoise = NULL,memb=u,mobileCenters=centers, fixedCenters=NULL, dist2clusters=dist2cent, withinss = withinss, size=size, functional=functional)
   class(res)<-"vegclust"
	return(res)
}

