clustcentroid <-
function(x,y) {
	if(is.null(dim(y))) {
	   u = as.memb(y)
	} else {
		u = as.matrix(y)
	}
   s = t(as.matrix(x))%*%u
   centers = sweep(t(s),1,colSums(u),"/")
   centers = as.data.frame(centers)   
   names(centers) = names(x)
   if(is.null(dim(y))) row.names(centers) = levels(as.factor(y))
   else row.names(centers) = names(y)
   
	return(centers)
}

