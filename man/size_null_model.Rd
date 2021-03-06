% Generated by roxygen2 (4.0.2): do not edit by hand
\name{size_null_model}
\alias{size_null_model}
\title{Size Ratio}
\usage{
size_null_model(speciesData, algo = "size_uniform", metric = "var_ratio",
  nReps = 1000, rowNames = TRUE, saveSeed = FALSE, algoOpts = list(),
  metricOpts = list())
}
\arguments{
\item{speciesData}{a dataframe <put some guidelines in here>}

\item{algo}{the algorithm to use, must be "size_uniform", "size_uniform_user", "size_source_pool", "size_gamma"}

\item{metric}{the metric used to caluclate the null model: choices are "min_diff", "min_ratio", "var_diff", "var_ratio"; default is Var.Ratio}

\item{nReps}{the number of replicates to run the null model.}

\item{rowNames}{Does your dataframe have row names? If yes, they are stripped, otherwise FALSE for data that has no row names}

\item{saveSeed}{TRUE or FALSE.  If TRUE the current seed is saved so the simulation can be repeated}

\item{algoOpts}{a list containing all the options for the specific algorithm you want to use.  Must match the algorithm given in the `algo` argument}

\item{metricOpts}{a list containing all the options for the specific metric you want to use.  Must match the metric given in the `metric` argument}
}
\description{
Create a size Ratio null model
}
\examples{
\dontrun{
## Run the null model
rodentMod <- size_null_model(dataRodents)
## Summary and plot info
summary(rodentMod)
plot(rodentMod,type="hist")
plot(rodentMod,type="size")

##  Uniform Size model with user inputs
rodentMod2 <- size_null_model(dataRodents,algo="size_uniform_user",
algoOpts = list(userLow = 3,userHigh=15))
summary(rodentMod2)
plot(rodentMod2,type="hist")
plot(rodentMod2,type="size")

### Source pool model

rodentMod_sp <- size_null_model(dataRodents,algo="size_source_pool",
algoOpts = list(sourcePool = runif(dim(dataRodents)[1],1,15)))

summary(rodentMod_sp)
plot(rodentMod_sp,type="hist")
plot(rodentMod_sp,type="size")

}
}

