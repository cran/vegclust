\encoding{UTF-8}
\name{plot.CAP}
\alias{plot.CAP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots cumulative abundance profiles
}
\description{
Create plots used to inspect one or more cumulative abundance profiles.
}
\usage{
\method{plot}{CAP}(x, sizes=NULL, species=NULL, plots=NULL, switchAxes=FALSE, add=FALSE, 
                   xlab="", ylab="",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object returned from function \code{\link{CAP}}.}
  \item{sizes}{A vector containing the (maximum) size values associated to each size class. If \code{NULL} the y-axis will be defined using the size class order in \code{x}.}
  \item{species}{A vector of strings indicating the species whose profile is to be drawn. If \code{NULL} all species are plotted.}
  \item{plots}{A vector of strings indicating the plot records whose profile is to be drawn. If \code{NULL} all plot records are plotted.}
  \item{switchAxes}{A flag indicating whether ordinate and abscissa axes should be interchanged.}
  \item{add}{A flag indicating whether profiles should be drawn on top of current drawing area. If \code{add=FALSE} a new plot is created.}
  \item{xlab}{String label for the x axis.}
  \item{ylab}{String label for the y axis.}
  \item{...}{Additional plotting parameters.}
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (in preparation) Dissimilarity measurements and the sie structure of ecological communities.

}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}
\seealso{
\code{\link{CAP}}
}
\examples{
## Load stratified data
data(medreg)

## Check that 'medreg' has correct class
class(medreg)

## Create cumulative abundance profile (CAP) for each plot
medreg.CAP = CAP(medreg)

## Draw the profile corresponding to the third plot
plot(medreg.CAP, plots="3")

## Look at the profile itself
medreg.CAP[[3]]
}
