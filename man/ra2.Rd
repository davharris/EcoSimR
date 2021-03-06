% Generated by roxygen2 (4.0.2): do not edit by hand
\name{ra2}
\alias{ra2}
\title{RA2 Niche Overlap Randomization Algorithm}
\usage{
ra2(speciesData = matrix(rpois(80, 1), nrow = 10))
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
replacing all non-zero elements with a random uniform [0,1] value.
}
\details{
The resource utilization matrix (rows = species, columns = discrete
resource categories) may include zeroes, but no negative numbers or missing
values. Relative resource within a species is first calculated, so the rows
need not sum to 1.0.
}
\note{
This algorithm retains the number and position of zero states in the
original matrix. However, all non-zero values are again replaced by a random
[0,1] value, which tends to inflate niche breadths of the simulated
assemblage. Although the results are not as severe as for RA1, this
algorithm is still prone to Type I errors, and is not recommended for niche
overlap analysis.
}
\examples{
ranUtil <- ra2(speciesData=matrix(rpois(40,0.5),nrow=8))
}
\references{
Kobayashi, S. 1991. Interspecific relations in forest floor
coleopteran assemblages: niche overlap and guild structure. Researches
in Population Ecology 33: 345-360.

Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural
assemblages of desert lizards and tropical fishes. Ecological Monographs
60: 27-55.
}

