clustmedoid <-
function(x, y, m=1) {
	if(is.null(dim(y))) {
	  u = as.memb(y)
	} else {
		u = as.matrix(y)^m
	}
  c = ncol(u)
  med = numeric(c)
  if(!inherits(x,"dist")) {
    d = as.matrix(dist(x))    
  }
	for(k in 1:c) {
	  med[k] = which.min((u[,k]^m)%*%d)
	}   
  names(med) = row.names(as.matrix(d))[med]
	return(med)
}

