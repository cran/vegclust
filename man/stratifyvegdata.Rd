\encoding{UTF-8}
\name{stratifyvegdata}
\alias{stratifyvegdata}
\alias{as.stratifiedvegdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reshape vegetation data into strata
}
\description{
Function \code{stratifyvegdata} reshapes individual plant abundance values into species abundance values per stratum. Function \code{as.stratifiedvegdata} checks if the input list has appropriate properties and turns it into an object of class '\code{stratifiedvegdata}'.
}
\usage{
stratifyvegdata(x, heights, plotColumn="plot", speciesColumn = "species", abundanceColumn="abundance", heightColumn = "height", counts=FALSE, mergeSpecies=FALSE)
as.stratifiedvegdata(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data frame containing individual plant data. Individuals are in rows, while measurements are in columns.}
  \item{heights}{A numerical vector containing the minimum height for each stratum in ascending order.}
  \item{plotColumn}{The name of the column in \code{x} that contains plot identifiers.}
  \item{speciesColumn}{The name of the column in \code{x} that contains species names.}
  \item{abundanceColumn}{The name of the column in \code{x} that contains abundance values.}
  \item{heightColumn}{The name of the column in \code{x} that contains height values.}
  \item{counts}{A flag to indicate that the output should be individual counts instead of added abundance values.}
  \item{mergeSpecies}{A flag to indicate that species identity should be ignored. This leads to analyzing the structure of biomass disregarding species identity.}
  \item{X}{A list with as many elements as plot records. Each element should be of class 'matrix' or 'data.frame' with species in rows and strata in columns. Furthermore, the number of rows (species) and columns (strata) should be the same for all elements.}
}
\details{
For each tree individual (row) in \code{x}, \code{stratifyvegdata} assigns it to the highest stratum whose height is lower or equal than the tree height. The corresponding abundance value is added to the abundance of the corresponding species at the stratum.
}
\value{
Both functions return an object of class '\code{stratifiedvegdata}', which is a list of data frames, one for each plot record. Each element (data frame) has as many rows as species and as many columns as strata (i.e., as many as elements in \code{height}). Columns are named starting with 'S' and continuing with the stratum number. If \code{mergeSpecies=TRUE} then all data frames have a single row (whose name is \code{"all"}).
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (in preparation) Dissimilarity assessments and the vertical structure of plant communities.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia.
}
\seealso{
\code{\link{reshape}}, \code{\link{CAP}}
}
\examples{
#Load tree data
data(treedata)

#Inspect tree data
head(treedata)

#Define stratum thresholds (4 strata)
heights = 0.5+(0:3)

#Stratify tree data
X = stratifyvegdata(treedata,heights=heights, plotColumn="plotID", speciesColumn="species", heightColumn="height", counts=TRUE)

#Inspect first plot
X[[1]]

}