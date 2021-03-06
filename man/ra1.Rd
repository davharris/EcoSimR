% Generated by roxygen2 (4.0.2): do not edit by hand
\name{ra1}
\alias{ra1}
\title{RA1 Niche Overlap Randomization Algorithm}
\usage{
ra1(speciesData = matrix(rpois(80, 1), nrow = 10))
}
\arguments{
\item{m}{a resource utilization matrix (rows = species, columns = discrete
resource states) filled with non-negative real numbers.}
}
\value{
Returns a random utilization matrix with the same dimensions as the
input matrix.
}
\description{
Randomizes a numeric utilization matrix speciesData by
replacing all elements with a random uniform [0,1] value.
}
\details{
The resource utilization matrix (rows = species, columns = discrete
resource categories) may include zeroes, but no negative numbers or missing
values. Relative resource within a species is first calculated, so the rows
need not sum to 1.0.
}
\note{
Because all matrix elements, including zeroes, are replaced with a
random uniform distribution, the null expectation is based on an assemblage
of generalist species with maximum niche breadth. This algorithm retains
neither the niche breadth of the individuals species nor the placement of 0
values (= unutilized resource states) in the matrix. These assumptions are
unrealistic, and a random matrix with zeroes will generate significantly low
niche overlap values with this metric. It is not recommended for niche
overlap analysis.
}
\examples{
ranUtil <- ra1(speciesData=matrix(rpois(40,0.5),nrow=8))
}
\references{
Kobayashi, S. 1991. Interspecific relations in forest floor
coleopteran assemblages: niche overlap and guild structure. Researches
in Population Ecology 33: 345-360.
}

