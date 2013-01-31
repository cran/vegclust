\encoding{UTF-8}
\name{vegdiststruct}
\alias{vegdiststruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Structural and compositional dissimilarity
}
\description{
Function to calculate the dissimilarity between plot records taking into account both vertical structure and composition.
}
\usage{
vegdiststruct(x, type="profile", method="bray", transform=NULL, stratumWeights=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A stratified vegetation data set (see function \code{\link{stratifyvegdata}}), or a set of cummulative abundance profiles (see function \code{\link{CAP}}).}
  \item{type}{Whether dissimilarities between pairs of sites should be calculated from differences in cumulative abundance profiles (\code{"profile"}), in total species abundance (\code{"abundance"}) or in volumes of cumulative abundance profiles (\code{"volume"}).}
  \item{method}{The dissimilarity index to calculate: (1) \code{method="bray"} for percentage difference (Bray-Curtis); (2) \code{method="ruzicka"} for Ruzicka index (a generalization of Jaccard); (3) \code{method="kulczynski"} for the Kulczynski dissimilarity index; and (4) \code{method="ochiai"} for the complement of a generalization of Ochiai index of similarity.}
  \item{transform}{A function or the name of a function to be applied to each cumulative abundance value.}
  \item{stratumWeights}{A numerical vector containing the weight for each stratum. If \code{NULL}, then all strata have the same weight.}
}
\details{
The following overview gives first the quantitative version, where y[ij] y[ik] refer to the quantity on column i and sites (rows) j and k of matrix Y. The four indices allowed are: 
\describe{
  \item{\code{method="bray"}:}{d[jk] = 1-2*(sum_i min(y[ij],y[ik]))/(yj_sum+yk_sum)}
  \item{\code{method="ruzicka"}:}{d[jk] = 1-(sum_i min(y[ij],y[ik]))/(sum_i max(y[ij],y[ik]))}
  \item{\code{method="kulczynski"}:}{d[jk] = 1-0.5*((sum_i min(y[ij],y[ik]))/yj_sum) + (sum_i min(y[ij],y[ik]))/yk_sum)}
  \item{\code{method="ochiai"}:}{d[jk] = 1-(sum_i min(y[ij],y[ik]))/sqrt(yj_sum*yk_sum)}
}
}
\value{
Returns an object of class '\code{\link{dist}}'.
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (in preparation) Dissimilarity assessments and the vertical structure of plant communities.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia.
}
\seealso{
\code{\link{stratifyvegdata}}, \code{\link[vegan]{vegdist}}
}
\examples{
#Load stratified data
data(medreg)

#Check that 'medreg' has correct class
class(medreg)

#Create cumulative abundance profile (CAP) for each plot
medreg.CAP = CAP(medreg)

#Create dissimilarity (percentage difference) matrix using profiles
medreg.D = vegdiststruct(medreg, method="bray")

#Create dissimilarity (percentage difference) matrix using abundances
medreg.D2 = vegdiststruct(medreg, method="bray", type="abundance")

#Calculate correlation
cor(as.vector(medreg.D), as.vector(medreg.D2))
}
