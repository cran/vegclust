clustcentroid <-
function(x,y, m=1) {
	if(is.null(dim(y))) {
	  u = as.memb(y)
	} else {
		u = as.matrix(y)^m
	}
   s = t(as.matrix(x))%*%(u)
   centers = sweep(t(s),1,colSums(u),"/")
   centers = as.data.frame(centers)   
   names(centers) = names(x)
   if(is.null(dim(y))) row.names(centers) = levels(as.factor(y))
   else row.names(centers) = names(y)
   
	return(centers)
}

