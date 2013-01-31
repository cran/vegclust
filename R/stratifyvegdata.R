stratifyvegdata<-function(x,heights, plotColumn="plot", speciesColumn = "species", abundanceColumn="abundance", heightColumn = "height", counts=FALSE, mergeSpecies=FALSE) {
  treeData = as.data.frame(x)
  plotColumnId = which(names(treeData)==plotColumn)
  abundanceColumnId = which(names(treeData)==abundanceColumn)
  heightColumnId = which(names(treeData)==heightColumn)
  speciesColumnId = which(names(treeData)==speciesColumn)
  if(mergeSpecies) treeData[,speciesColumnId] = "allspecies"
  spnames =unique(treeData[,speciesColumnId])
 
  stratify<-function(treeDataPlot, heights, spnames=NULL, speciesColumnId, abundanceColumnId, heightColumnId, counts=FALSE) {
    if(is.null(spnames)) spnames = unique(treeData[,speciesColumnId])
    nsp = length(spnames)
    nstrata = length(heights)
    m = data.frame(matrix(0,nrow=nsp, ncol=nstrata))
    row.names(m) = spnames
    if(!is.null(names(heights))) names(m) = names(heights)
    else names(m) = paste("S",1:nstrata, sep="")
    for(i in 1:nrow(treeDataPlot)) {
      isp = which(spnames==treeDataPlot[i,speciesColumnId])
      sel = (heights <= treeDataPlot[i, heightColumnId])
      if(sum(sel,na.rm=T)>0) {
        istratum = max(which(sel))
        if(!counts) m[isp,istratum] = m[isp,istratum]+treeDataPlot[i,abundanceColumnId]
        else m[isp,istratum] = m[isp,istratum]+1
      }
    }
    return(m)
  }
  
  X = lapply(split(treeData,treeData[,plotColumnId]),
         FUN=stratify,heights=heights, 
         spnames=spnames, speciesColumnId = speciesColumnId, abundanceColumnId =abundanceColumnId, heightColumnId=heightColumnId, counts=counts)
  class(X)<-c("list","stratifiedvegdata")
  return(X)
}