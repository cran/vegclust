\name{treedata}
\alias{treedata}
\docType{data}
\title{Synthetic vegetation data set with tree data}
\description{
  A synthetic data set used to illustrate the stratification of data originally collected on an individual basis.
}
\usage{data(treedata)}
\format{
  A data frame where each row corresponds to a different tree. Columns are plot code, species identity, tree height and cover value.
}
\seealso{\code{\link{stratifyvegdata}}}
\examples{
#Load tree data
data(treedata)

#Inspect tree data
head(treedata)

#Define stratum thresholds (4 strata)
heights = 0.5+(0:3)

#Stratify tree data
X = stratifyvegdata(treedata,heights=heights, plotColumn="plotID", speciesColumn="species", heightColumn="height", counts=TRUE)

#Inspect first plot
X[[1]]
}
\keyword{datasets}
