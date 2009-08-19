\encoding{UTF-8}
\name{relate.levels}
\alias{relate.levels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relates two clustering level results.
}
\description{
Analyzes how lower level clusters are assigned into upper level ones. The analysis is made for several number of clusters.
}
\usage{
relate.levels(lower, upper, verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower}{ A list of objects of type \code{\link{vegclust}} that represent classifications at a fine level of resolution.}
  \item{upper}{ A list of objects of type \code{\link{vegclust}} that represent classifications at an broad level of resolution.}
  \item{verbose}{ A flag used to ask for extra screen output.}
  \item{...}{
    Additional parameters for function \code{\link{defuzzify}}.
  }
}
\details{
For each pair of \code{vegclust} objects in \code{upper} and \code{lower}, the function makes calls to \code{\link{vegclass}} in order to assign the centroids of the object of \code{lower} into classes of the object of \code{upper}. The assignment result is then defuzzified (by calling function \code{\link{defuzzify}}) and several quantities are calculated (see 'value' section).
}
\value{
 A list with several data frames (see below). In each of them, the rows are items of \code{upper} and columns are items of \code{lower}. The names of rows and columns are the number of clusters of each \code{\link{vegclust}} object.
 \item{inter}{The number of low level clusters that are not clearly assigned into any cluster (with the options supplied to \code{\link{defuzzify}}). }
 \item{noise}{The number of low level clusters that are assigned to the Noise class (for \code{upper} objects using Noise clustering). }
 \item{miss}{The number of upper level clusters (mobile or fixed) that do not have any member among the lower level clusters. }
 \item{missmob}{The number of upper level mobile clusters that do not have any member among the lower level clusters. }
 \item{missfix}{The number of upper level fixed clusters that do not have any member among the lower level clusters. }
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}

\seealso{
\code{\link{vegclust}}, \code{\link{vegclass}}, \code{\link{defuzzify}}
}
\examples{
# Loads data  
data(wetland)
  
# This equals the chord transformation (see also 'normalize' option in \code{\link{decostand}} from the vegan package)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, sqrt(rowSums(as.matrix(wetland)^2)), "/"))

# Create noise clustering from hierarchical clustering at different number of cluster
wetland.hc = hclust(dist(wetland.chord),method="ward") 
wetland.nc1 = hier.vegclust(wetland.chord, wetland.hc, cmin=2, cmax=6, m = 1.2, dnoise=0.75, method="NC")
wetland.nc2 = hier.vegclust(wetland.chord, wetland.hc, cmin=2, cmax=4, m = 1.2, dnoise=0.85, method="NC")

#Studies the assignment of levels
relate.levels(wetland.nc1, wetland.nc2, method="cut")
}
