% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkSolution.R
\name{checkSolution}
\alias{checkSolution}
\title{Checks whether the solution generated by the ODE solver contains negative values}
\usage{
checkSolution(soln, tol = -0.1)
}
\arguments{
\item{soln}{Matrix from ode solver out$solution}

\item{tol}{tolerance}
}
\description{
Checks whether the solution generated by the ODE solver contains negative values
}
