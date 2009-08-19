\encoding{UTF-8}
\name{clustcentroid}
\alias{clustcentroid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Centroids of a crisp classification
}
\description{
Calculates the centroid (multivariate average) coordinates of a classification represented in the form of a cluster vector or membership matrix.
}
\usage{
clustcentroid(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Community data, a site-by-species data frame.
   }
  \item{y}{
	A vector indicating the cluster that each object in \code{x} belongs to. Alternatively, a fuzzy/hard site-by-group matrix of membership values.
   }
}
\value{
 A group-by-species data frame containing species average abundance values.
}
\note{
  In order to assign new plot record data into a predefined set of classes, one should use functions \code{\link{as.vegclust}} and \code{\link{vegclass}} instead.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}

\seealso{
\code{\link{as.vegclust}}, \code{\link{vegclass}}, \code{\link{vegclust}}, \code{\link{kmeans}}
}
\examples{
# Loads stats
library(stats)
  
# Loads data
data(wetland)
  
# This equals the chord transformation (see also 'normalize' option in \code{\link{decostand}} from the vegan package)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, sqrt(rowSums(as.matrix(wetland)^2)), "/"))

# Performs a K-means clustering
wetland.km = kmeans(wetland.chord, centers=3, nstart=10)

# Gets cluster variance of KM (i.e. hard) clusters
clustcentroid(wetland.chord, y=wetland.km$cluster)

}
