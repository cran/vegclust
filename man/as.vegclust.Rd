\encoding{UTF-8}
\name{as.vegclust}
\alias{as.vegclust}
\title{
Turns into vegclust objects
}
\description{
Attempts to turn its arguments into a \code{\link{vegclust}} object
}
\usage{
as.vegclust(x, cluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A site-by-species data matrix (raw mode), or a site-by-site distance matrix (distance mode).
   }
  \item{cluster}{
	A vector indicating the hard membership of each object in \code{x} to a set of groups.
   }
}
\details{
This function is used to generate \code{\link{vegclust}} objects which can then be used in \code{\link{vegclass}} to classify new data. The input classification is hard (i.e. yes/no membership), cluster centers are calculated as multivariate means, and the method for assigning new data is assumed to be k-means (\code{"KM"}), i.e. plots will be assigned to the nearest cluster center. If community data is given as site-by-species data matrix the cluster centroids are added as \code{mobileCenters} in the \code{\link{vegclust}} object. Centroids cannot be computed if community data is given as a site-by-site dissimilarity matrix. 
}
\value{
An object of class \code{\link{vegclust}}.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia.
}

\seealso{
\code{\link{vegclust}}, \code{\link{vegclass}}
}
\examples{
# Loads data  
data(wetland)

# This equals the chord transformation (see also 'normalize' option in \code{\link{decostand}} from the vegan package)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, sqrt(rowSums(as.matrix(wetland)^2)), "/"))

# Splits wetland data into two matrices of 30x27 and 11x22
wetland.30 = wetland.chord[1:30,]
wetland.30 = wetland.30[,colSums(wetland.30)>0]
dim(wetland.30)
wetland.11 = wetland.chord[31:41,]
wetland.11 = wetland.11[,colSums(wetland.11)>0] 
dim(wetland.11)

# Performs a K-means clustering of the data set with 30 sites
wetland.km = kmeans(wetland.30, centers=3, nstart=10)

# Transforms the 'external' classification of 30 sites into a 'vegclust' object
wetland.30.vc = as.vegclust(wetland.30, wetland.km$cluster)

# Assigns the second set of sites according to the (k-means) membership rule 
# That is, sites are assigned to the cluster whose cluster centroids is nearest.
wetland.11.km = vegclass(wetland.30.vc, wetland.11)

# A similar 'vegclust' object is obtained when using the distance mode...
wetland.d.vc = as.vegclust(dist(wetland.30), wetland.km$cluster)

# which can be also used to produce the assignment of the second set of objects
wetland.d.11 = as.data.frame(as.matrix(dist(wetland.chord)))[31:41,1:30]
wetland.d.11.km = vegclass(wetland.d.vc,wetland.d.11)
}
