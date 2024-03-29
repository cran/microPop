% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/productionFuncDefault.R
\name{productionFuncDefault}
\alias{productionFuncDefault}
\alias{productionFunc}
\title{Production Function}
\usage{
productionFuncDefault(
  strainName,
  groupName,
  pathName,
  varName,
  all.substrates,
  keyResName,
  stoichiom,
  products,
  bio.products,
  uptake,
  growthRate,
  yield,
  parms,
  water
)
}
\arguments{
\item{strainName}{Name of the strain that is being looped through in the ODE solver}

\item{groupName}{Name of microbial group that is being looped through in the ODE solver}

\item{pathName}{Name of metabolic path (e.g. path1) that is being looped through in the ODE solver}

\item{varName}{(string). Calculate production of this variable}

\item{all.substrates}{Vector of strings giving the names of the all the substrates used on this pathway}

\item{keyResName}{(string). Name of the key resource on this pathway}

\item{stoichiom}{Named vector (names are resourceNames) giving the mass of each resource in the stoichiometry i.e. molar mass of resource multiplied by the number of moles in the stoichiometry}

\item{products}{Vector of strings giving the names of the all the metabolic products created on this pathway}

\item{bio.products}{Vector of strings giving the names of the all the microbial products created on this pathway}

\item{uptake}{Vector with names given by resourceNames which given mass uptake of each resource per unit time}

\item{growthRate}{(scalar) microbial growth rate (mass per unit time) on the given pathway}

\item{yield}{Named vector (names are resourceNames) giving the mass yield of biomass on each resource (mass microbe/mass resource)}

\item{parms}{List containing all system parameters}

\item{water}{Name of resource with Rtype 'Sw' - i.e resource could be called 'water' or 'H2O' etc}
}
\value{
(scalar) production rate of given resource (units are resource mass/time)
}
\description{
Production rate of resource (units are resource mass/time)
}
