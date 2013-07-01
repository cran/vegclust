\encoding{UTF-8}
\name{vegdiststruct}
\alias{vegdiststruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Structural and compositional dissimilarity
}
\description{
Function to calculate the dissimilarity between plant communities taking into account both their composition and the size of organisms.
}
\usage{
vegdiststruct(x, type="profile", method="bray", transform=NULL, classWidths=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A stratified vegetation data set (see function \code{\link{stratifyvegdata}}), or a set of cummulative abundance profiles (see function \code{\link{CAP}}).}
  \item{type}{Whether dissimilarities between pairs of sites should be calculated from differences in cumulative abundance profiles (\code{"profile"}), in total species abundance (\code{"abundance"}) or in volumes of cumulative abundance profiles (\code{"volume"}).}
  \item{method}{The dissimilarity coefficient to calculate (see details).}
  \item{transform}{A function or the name of a function to be applied to each cumulative abundance value.}
  \item{classWidths}{A numerical vector containing the width for each size class (or stratum). If \code{NULL}, then the function assumes classes of equal width.}
}
\details{
The six different coefficients available are described in De Caceres et al. (in prep.): (1) \code{method="bray"} for percentage difference (alias Bray-Curtis dissimilarity); (2) \code{method="ruzicka"} for Ruzicka index (a generalization of Jaccard); (3) \code{method="kulczynski"} for the Kulczynski dissimilarity index; (4) \code{method="ochiai"} for the complement of a quantitative generalization of Ochiai index of similarity; (5) \code{method="canberra"} for the Canberra index (Adkins form); (6) \code{method="relman"} for the relativized Manhattan coefficient (Whittaker's index of association).
}
\value{
Returns an object of class '\code{\link{dist}}'.
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (in preparation) Dissimilarity measurements and the sie structure of ecological communities.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia.
}
\seealso{
\code{\link{stratifyvegdata}}, \code{\link[vegan]{vegdist}}
}
\examples{
## Load stratified data
data(medreg)

## Check that 'medreg' has correct class
class(medreg)

## Create cumulative abundance profile (CAP) for each plot
medreg.CAP = CAP(medreg)

## Create dissimilarity (percentage difference) matrix using profiles
medreg.D = vegdiststruct(medreg, method="bray")

## Create dissimilarity (percentage difference) matrix using abundances
medreg.D2 = vegdiststruct(medreg, method="bray", type="abundance")

## Calculate correlation
cor(as.vector(medreg.D), as.vector(medreg.D2))

}
