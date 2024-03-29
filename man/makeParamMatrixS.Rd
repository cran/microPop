% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeParamMatrixS.R
\name{makeParamMatrixS}
\alias{makeParamMatrixS}
\title{Gets parameter values for parameters halfSat, yield and maxGrowthRate
from the MFGs and puts into a matrix
also assigns strain traits
since traits are assigned to all strains at once for one param they are
stored in Mat[strain,res,path], this is then rearranged to make a matrix
for each strain (mat[path,res]).}
\usage{
makeParamMatrixS(
  resNames,
  microbeNames,
  parameterName,
  numPaths,
  numStrains,
  strainOptions,
  oneStrainRandomParams
)
}
\arguments{
\item{resNames}{Vector of strings which contains the names of the resources in the system}

\item{microbeNames}{Vector of strings which contains the names of the microbial groups in the system e.g. c('Bacteroides','Acetogens')}

\item{parameterName}{Name of parameter}

\item{numPaths}{Named vector. Number of paths for each microbial group}

\item{numStrains}{Integer or named vector of integers. Number of strains per group}

\item{strainOptions}{List of strain options}

\item{oneStrainRandomParams}{Logical. TRUE for randomized params even if there is only one strain.}
}
\value{
A list called parameterName which contains matrices for all strains (in all groups) with paths on rows and resources on columns
}
\description{
Gets parameter values for parameters halfSat, yield and maxGrowthRate
from the MFGs and puts into a matrix
also assigns strain traits
since traits are assigned to all strains at once for one param they are
stored in Mat[strain,res,path], this is then rearranged to make a matrix
for each strain (mat[path,res]).
}
\keyword{internal}
