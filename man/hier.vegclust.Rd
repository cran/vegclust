\encoding{UTF-8}
\name{hier.vegclust}
\alias{hier.vegclust}
\alias{random.vegclust}
\title{
Clustering with several number of clusters
}
\description{
Performs several runs of function 'vegclust' on a community data matrix using different number of clusters
}
\usage{
hier.vegclust(x, hclust, cmin=2,cmax=20, verbose=TRUE, ...)
random.vegclust(x, cmin=2, cmax=20, nstart=10, verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	Community data table. A site (rows) by species (columns) matrix or data frame.
  }
  \item{hclust}{
	A hierarchical clustering represented in an object of type \code{\link{hclust}}.
  }
  \item{cmin}{
	Number of minimum mobile clusters.
   }
  \item{cmax}{
	Number of maximum mobile clusters.
   }
  \item{nstart}{
	A number indicating how many random trials should be performed for each number of groups
  }
  \item{verbose}{
 	Flag used to print which number of clusters is currently running.
  }
  \item{...}{
	Additional parameters for function \code{\link{vegclust}}.
  }
}
\details{
Function \code{hier.vegclust} takes starting cluster configurations from cuts of a dendrogram given by object \code{hclust}. Function \code{random.vegclust} chooses random objects as cluster centroids and for each number of clusters performs \code{nstart} trials.
}
\value{
Returns an object of type 'mvegclust' (multiple vegclust), which contains a list vector with objects of type \code{\link{vegclust}}.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}


\seealso{
\code{\link{vegclust}},\code{\link{vegclass}},\code{\link{defuzzify}},\code{\link{hclust}}
}
\examples{

# Loads data  
data(wetland)
  
# This equals the chord transformation (see also 'normalize' option in \code{\link{decostand}} from the vegan package)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, sqrt(rowSums(as.matrix(wetland)^2)), "/"))

# Create noise clustering from hierarchical clustering at different number of clusters
wetland.hc = hclust(dist(wetland.chord),method="ward") 
wetland.nc1 = hier.vegclust(wetland.chord, wetland.hc, cmin=2, cmax=5, m = 1.2, dnoise=0.75, method="NC")

# Create noise clustering from random seeds at different levels
wetland.nc2 = random.vegclust(wetland.chord, cmin=2, cmax=5, nstart=10, m = 1.2, dnoise=0.75, method="NC")
}
