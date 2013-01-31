\encoding{UTF-8}
\name{CAP}
\alias{CAP}
\alias{CAP2matrix}
\alias{CAPmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cummulative abundance profile (CAP)
}
\description{
Functions to calculate cummulative abundance profiles (CAPs), to build matrices from them, and to average profiles.
}
\usage{
CAP(x, transform=NULL)
CAP2matrix(CAP, type="profile",stratumWeights=NULL)
CAPmeans(CAP, factor=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A stratified vegetation data set (see function \code{\link{stratifyvegdata}}).}
  \item{transform}{A function or the name of a function to be applied to each cumulative abundance value.}
  \item{CAP}{An object of class '\code{CAP}'.}
  \item{type}{The type of information that the resulting matrix should contain. Either \code{"profile"}, \code{"abundance"} or \code{"volume"}.}
  \item{stratumWeights}{A numerical vector containing the weight for each stratum. If \code{NULL}, then all strata have the same weight.}
  \item{factor}{A factor to calculate average profiles per each level.}
}
\details{
Function \code{CAP} replaces the abundance value of a stratum by the sum of abundances in this and higher strata. Thus, upper strata contain smaller abundance values than lower strata, creating a cumulative abundance profile. Function \code{CAP2matrix} taks an object of class '\code{CAP}' and returns a data matrix, where values differ depending on parameter \code{type}: (1) \code{type="profile"} simply reshapes the '\code{CAP}' object (a list) into a matrix with as many rows as plot records and where columns are organized in blocks (there are as many blocks as species and each block has as many columns as strata); (2) \code{type="abundance"} returns a plot-by-species matrix where each value is the total abundance of the species in the plot (i.e. the CAP value at the ground level); (3) \code{type="volume"} returns a plot-by-species matrix where each value is the sum of CAP values across strata (a measure of the "volume" occupied by the species in the plot). When provided, \code{stratumWeights} are used to weight strata of the cumulative abundance profiles (for (1) and (3) only). Function \code{CAPmeans} calculates the average abundance profile for a set of plot records. If supplied, \code{factor} can be used to speficy groups of samples for which average profiles are to be calculated. 
}
\value{
Function \code{CAP} returns an object of class '\code{CAP}', similar to objects of class '\code{stratifiedvegdata}' but where abundance values of upper strata have beed added to those of lower strata. Function \code{CAP2matrix} returns a matrix with species as rows (columns depend on the value of \code{type}).
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (in preparation) Dissimilarity assessments and the vertical structure of plant communities.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia.
}
\seealso{
  \code{\link{stratifyvegdata}}, \code{\link{plot.CAP}}, \code{\link{vegdiststruct}}
}
\examples{
#Load stratified data
data(medreg)

#Check that 'medreg' has correct class
class(medreg)

#Look at the data for the third plot
medreg[[3]]

#Create cumulative abundance profile (CAP) for each plot
medreg.CAP = CAP(medreg)

#Look at the profile of the third plot
medreg.CAP[[3]]

#Create matrix with species abundances
medreg.X = CAP2matrix(medreg.CAP, type="abundance")
head(medreg.X)

#Generate average profile
average.CAP = CAPmeans(medreg.CAP)

#Plot average profile
plot(average.CAP)

}
