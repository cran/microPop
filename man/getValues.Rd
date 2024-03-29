% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getValues.R
\name{getValues}
\alias{getValues}
\title{get system quantity (e.g. startValue, inflowRate, washOut) for all state variables
(convention is that microbes are before resources)}
\usage{
getValues(
  sysInfoMicrobes,
  sysInfoRes,
  stateVarNames,
  quantity,
  strainNames,
  microbeNames,
  resourceNames,
  numStrains
)
}
\arguments{
\item{sysInfoMicrobes}{sys info dataframe for microbes}

\item{sysInfoRes}{sys info dataframe for resources}

\item{stateVarNames}{Vector of names of all the state variables}

\item{quantity}{String. Name of quantity to get value for e.g. 'startValue'}

\item{strainNames}{Vector of strings of strain names}

\item{microbeNames}{Vector of strings which contains the names of the microbial groups in the system e.g. c('Bacteroides','Acetogens')}

\item{resourceNames}{Vector of strings which contains the names of the resources in the system}

\item{numStrains}{Integer. Number of strains per group}
}
\description{
get system quantity (e.g. startValue, inflowRate, washOut) for all state variables
(convention is that microbes are before resources)
}
