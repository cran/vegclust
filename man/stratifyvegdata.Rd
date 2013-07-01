\encoding{UTF-8}
\name{stratifyvegdata}
\alias{stratifyvegdata}
\alias{as.stratifiedvegdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reshapes community data from individual into stratified form
}
\description{
Function \code{stratifyvegdata} reshapes individual abundance values into species abundance values per size class (stratum). Function \code{as.stratifiedvegdata} checks if the input list has appropriate properties and turns it into an object of class '\code{stratifiedvegdata}'.
}
\usage{
stratifyvegdata(x, sizes, plotColumn="plot", speciesColumn = "species", 
                abundanceColumn="abundance", sizeColumn = "size", counts=FALSE, 
                mergeSpecies=FALSE)
as.stratifiedvegdata(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data frame containing individual plant data. Individuals are in rows, while measurements are in columns.}
  \item{sizes}{A numerical vector containing the minimum size for each size class (stratum) in ascending order.}
  \item{plotColumn}{The name of the column in \code{x} that contains plot identifiers.}
  \item{speciesColumn}{The name of the column in \code{x} that contains species names.}
  \item{abundanceColumn}{The name of the column in \code{x} that contains abundance values.}
  \item{sizeColumn}{The name of the column in \code{x} that contains size values.}
  \item{counts}{A flag to indicate that the output should be individual counts instead of added abundance values.}
  \item{mergeSpecies}{A flag to indicate that species identity should be ignored. This leads to analyzing the structure of biomass disregarding species identity.}
  \item{X}{A list with as many elements as plot records. Each element should be of class 'matrix' or 'data.frame' with species in rows and strata in columns. Furthermore, the number of rows (species) and columns (strata) should be the same for all elements.}
}
\details{
For each individual (row) in \code{x}, \code{stratifyvegdata} assigns it to the size class (stratum) containing its size. The corresponding abundance value (e.g. crown cover) of the individual is added to the abundance of the corresponding species at the size class (stratum).
}
\value{
Both functions return an object of class '\code{stratifiedvegdata}', which is a list of data frames, one for each plot record. Each element (data frame) has as many rows as species and as many columns as size classes (i.e., as many as elements in vector \code{sizes}). Columns are named starting with 'S' and continuing with the size class (stratum) number. If \code{mergeSpecies=TRUE} then all data frames have a single row (whose name is \code{"all"}).
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (in preparation) Dissimilarity measurements and the sie structure of ecological communities.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia.
}
\seealso{
\code{\link{reshape}}, \code{\link{CAP}}
}
\examples{
## Load tree data
data(treedata)

## Inspect tree data
head(treedata)

## Define stratum thresholds (4 strata)
heights = 0.5+(0:3)

## Stratify tree data
X = stratifyvegdata(treedata, sizes=heights, plotColumn="plotID",
                    speciesColumn="species", sizeColumn="height", counts=TRUE)

## Inspect the first plot record
X[[1]]

}