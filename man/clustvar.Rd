\encoding{UTF-8}
\name{clustvar}
\alias{clustvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cluster variance (beta diversity)
}
\description{
Computes the variation in community composition (i.e. beta diversity) found within the sites of a set of hard or fuzzy clusters.
}
\usage{
clustvar(x, cluster = NULL, defuzzify=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Community data. Either a site-by-species matrix or a site-by-site matrix of compositional distances between sites  (i.e., an object of class \code{\link{dist}}). Alternatively, this can be an object of class \code{\link{vegclust}} or \code{\link{vegclass}}, and in this case it is unnecessary to provide \code{cluster}.
   }
  \item{cluster}{
	A vector indicating the hard membership of each object in \code{x} to a set of groups.
   }
  \item{defuzzify}{
        A flag indicating whether fuzzy memberships should be defuzzified (see function \code{\link{defuzzify}}). Only applies to the case where an object of class \code{\link{vegclust}} or \code{\link{vegclass}} is supplied in \code{x}.
   }
  \item{...}{
    Additional parameters for function \code{\link{defuzzify}}.
  }
}
\details{
If an object of class \code{\link{vegclust}} or \code{\link{vegclass}} is provided as input, the function uses the information contained there (distances to clusters, memberships and exponent of fuzziness) in order to compute cluster variances. Both mobile and fixed clusters are considered. Membership matrices may be defuzzified if \code{defuzzify=TRUE}. Alternatively, if \code{x} is a data matrix (site by species or distances among sites) and \code{cluster} is \code{null}, the function assumes a single cluster of all points in \code{x}. When \code{cluster} is provided, the function computes cluster variance for each (hard) group and this computation implies setting the centroid of the group. In contrast, when an object of class \code{\link{vegclust}} or \code{\link{vegclass}} is provided the centroids do not need to be recomputed, and the distances to centroids are used directly.
}
\value{
 A double value (for one cluster) or a vector of values, one per each cluster.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}

\seealso{
\code{\link{vegclust}}, \code{\link{kmeans}},\code{\link{defuzzify}}
}
\examples{
# Loads data  
data(wetland)
  
# This equals the chord transformation (see also 'normalize' option in \code{\link{decostand}} from the vegan package)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, sqrt(rowSums(as.matrix(wetland)^2)), "/"))

# Create noise clustering with 3 clusters. Perform 10 starts from random seeds and keep the best solution
wetland.nc = vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, method="NC", nstart=10)

# Gets cluster variance of fuzzy clusters
clustvar(wetland.nc)

# Gets cluster variance of fuzzy clusters after defuzzification
clustvar(wetland.nc, defuzzify=TRUE)

# Similar to the previous, this gets cluster variance of defuzzified (i.e. hard) clusters
clustvar(wetland.chord, cluster=defuzzify(wetland.nc)$cluster)

# Gets cluster variance of K-means (i.e. hard) clusters
clustvar(wetland.chord, cluster=kmeans(wetland.chord, centers=3, nstart=10)$cluster)

}
