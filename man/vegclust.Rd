\encoding{UTF-8}
\name{vegclust}
\alias{vegclust}
\alias{vegclustdist}
\title{
Vegetation clustering methods
}
\description{
Performs hard or fuzzy clustering of a community data matrix
}
\usage{
vegclust(x, mobileCenters, fixedCenters = NULL, method="NC", m = 2, dnoise = NULL, eta = NULL, alpha=0.001, iter.max=100, nstart=1, maxminJ = 10, seeds=NULL, verbose=FALSE)
vegclustdist(x, mobileMemb, fixedMemb = NULL, method="NC", m = 2, dnoise = NULL, eta = NULL, alpha=0.001, iter.max=100, nstart=1, seeds=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	Community data. A site-by-species matrix or data frame (for \code{vegclust}) or a site-by-site dissimilarity matrix or \code{\link{dist}} object (for \code{vegclustdist}).
  }
  \item{mobileCenters}{
	A number, a vector of seeds, or coordinates for mobile clusters.
   }
  \item{fixedCenters}{
	A matrix or data frame with coordinates for fixed (non-mobile) clusters.
  }
  \item{mobileMemb}{
	A number, a vector of seeds, or starting memberships for mobile clusters.
   }
  \item{fixedMemb}{
	A matrix or data frame with starting memberships for fixed (non-mobile) clusters.
  }
  \item{method}{
	A clustering model. Current accepted models are: \code{"KM"} for kmeans (MacQueen 1967), \code{"FCM"} for fuzzy c-means (Bezdek 1981), \code{"NC"} for noise clustering (Dave and Krishnapuram 1997) and \code{"PCM"} for possibilistic c-means (Krishnapuram and Keller 1993).
  }
  \item{m}{
	The fuzziness exponent to be used (this is relevant for all models except for kmeans)
  }
  \item{dnoise}{
	The distance to the noise cluster, relevant for noise clustering (NC). 
  }
  \item{eta}{
	A vector of reference distances, relevant for possibilistic C-means (PCM). 
  }
  \item{alpha}{
	Threshold used to stop iterations. The maximum difference in the membership matrix of the current vs. the previous iteration will be compared to this value.
  }
  \item{iter.max}{
	The maximum number of iterations allowed.
  }
  \item{nstart}{
	If \code{mobileCenters} or \code{mobileMemb} is a number, how many random sets should be chosen?
  }
  \item{maxminJ}{
	When random starts are used, these will stop if at least \code{maxminJ} runs ended up in the same functional value.
  }
  \item{seeds}{
	If \code{mobileCenters} or \code{mobileMemb} is a number, a vector indicating which objects are potential initial centers. If \code{NULL} all objects are valid seeds.
  }
  \item{verbose}{
	Flag to print extra output.
  }

}
\details{
Functions \code{vegclust} and \code{vegclustdist} try to generalize the \code{\link{kmeans}} function in \code{stats} in three ways. Firstly, they allows different clustering models (three fuzzy and one hard). The reader should refer to the original publications to better understand the differences between models. Secondly, users can specify fixed clusters (that is, centroids that do not change their positions during iterations). Fixed clusters are intended to be used when some clusters were previously defined and new data has been collected. One may allow some of these new data points to form new clusters, while some other points will be assigned to the original clusters. In the case of models with cluster repulsion (such as KM, FCM or NC) the new (mobile) clusters are not allowed to 'push' the fixed ones. As a result, mobile clusters will occupy new regions of the reference space. Thirdly, \code{vegclustdist} implements the distance-based equivalent of \code{vegclust}. Note that all data frames or matrices used as input of \code{vegclust} should be defined on the same space of species (see \code{\link{conformveg}}). Unlike \code{\link{kmeans}}, which allows different specific algorithms, here updates of centroids are done after all objects have been reassigned (Forgy 1965). In order to obtain hard cluster definitions, users can apply the function \code{\link{defuzzify}} to the \code{vegclust} object.
}
\value{
Returns an object of type \code{vegclust} with the following items:
  \item{mode}{\code{raw} for function \code{vegclust} and \code{dist} for function \code{vegclustdist}.}
  \item{method}{The clustering model used}
  \item{m}{The fuzziness exponent used (\code{m=1} in case of kmeans)}
  \item{dnoise}{The distance to the noise cluster used for noise clustering (NC). This is set to \code{NULL} for other models.}
  \item{eta}{ The reference distance vector used for possibilistic c-means (PCM). This is set to \code{NULL} for other models.}
  \item{memb}{The fuzzy membership matrix. Columns starting with "M" indicate mobile clusters, whereas columns starting with "F" indicate fixed clusters.}
  \item{mobileCenters}{A data frame with the coordinates of the mobile centers (not available in \code{vegclustdist}).}
  \item{fixedCenters}{A data frame with coordinates for fixed (non-mobile) clusters (not available in \code{vegclustdist}).}
  \item{dist2clusters}{The matrix of object distances to cluster centers. Columns starting with "M" indicate mobile clusters, whereas columns starting with "F" indicate fixed clusters.}
  \item{withinss}{The within-cluster sum of squares for each cluster.}
  \item{size}{The number of objects belonging to each cluster. In case of fuzzy clusters the sum of memberships is given.}
  \item{functional}{The objective function value (the minimum value attained after all iterations).}
}
\references{
Forgy, E. W. (1965) Cluster analysis of multivariate data: efficiency vs interpretability of classifications. Biometrics 21, 768-769.

MacQueen, J. (1967) Some methods for classification and analysis of multivariate observations. In Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability, eds L. M. Le Cam and J. Neyman, 1, pp. 281-297. Berkeley, CA: University of California Press.

\enc{Davé}{Dave}, R. N. and R. Krishnapuram (1997) Robust clustering methods: a unified view. IEEE Transactions on Fuzzy Systems 5, 270-293.

Bezdek, J. C. (1981) Pattern recognition with fuzzy objective functions. Plenum Press, New York.

Krishnapuram, R. and J. M. Keller. (1993) A possibilistic approach to clustering. IEEE transactions on fuzzy systems 1, 98-110.

De \enc{Cáceres}{Caceres}, M., Font, X, Oliva, F. (2010) The management of numerical vegetation classifications with fuzzy clustering methods. Journal of Vegetation Science 21 (6): 1138-1151.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}


\seealso{
\code{\link{kmeans}}, \code{\link[e1071]{cmeans}},\code{\link{vegclass}},\code{\link{defuzzify}}
}
\examples{

# Loads data  
data(wetland)
  
# This equals the chord transformation (see also 'normalize' option in \code{\link{decostand}} from the vegan package)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, sqrt(rowSums(as.matrix(wetland)^2)), "/"))

# Create noise clustering with 3 clusters. Perform 10 starts from random seeds and keep the best solution
wetland.nc = vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, method="NC", nstart=10)

# Fuzzy membership matrix
wetland.nc$memb

# Cardinality of fuzzy clusters (i.e., the number of objects belonging to)
wetland.nc$size

# Obtains hard membership vector, with 'N' for objects that are unclassified
defuzzify(wetland.nc$memb)$cluster

# The same result is obtained with a matrix of chord distances
wetland.d = dist(wetland.chord)
wetland.d.nc = vegclustdist(wetland.d, mobileMemb=3, m = 1.2, dnoise=0.75, method="NC", nstart=10)

}
