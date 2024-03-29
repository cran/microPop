% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAllResources.R
\name{getAllResources}
\alias{getAllResources}
\title{Makes vector of unique resource names}
\usage{
getAllResources(microbeNames, gutModel = FALSE, myPars = NULL)
}
\arguments{
\item{microbeNames}{Vector of strings which contains the names of the microbial groups in the system e.g. c('Bacteroides','Acetogens')}

\item{gutModel}{Logical. TRUE if using with the microPopGut package}

\item{myPars}{list of extra parameters}
}
\value{
vector of resource names
}
\description{
Makes vector of unique resource names
}
