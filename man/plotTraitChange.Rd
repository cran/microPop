% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotTraitChange.R
\name{plotTraitChange}
\alias{plotTraitChange}
\title{plot changes in trait over time}
\usage{
plotTraitChange(
  out,
  trait.name,
  group.names,
  resource.name = NULL,
  path = NULL,
  xlabel = "Time (days)",
  saveFig = FALSE,
  figType = "eps",
  figName = "Traits"
)
}
\arguments{
\item{out}{Output from microPopModel()}

\item{trait.name}{can be 'halfSat','yield','maxGrowthRate' and 'pHtrait' or 'strainpHcorners'}

\item{group.names}{can be a vector of group names or just one string for one name}

\item{resource.name}{String}

\item{path}{String}

\item{xlabel}{String}

\item{saveFig}{Logical}

\item{figType}{String}

\item{figName}{String}
}
\description{
plot changes in trait over time
}
