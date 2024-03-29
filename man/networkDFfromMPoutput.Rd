% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networkDFfromMPoutput.R
\name{networkDFfromMPoutput}
\alias{networkDFfromMPoutput}
\title{networkDFfromMPoutput}
\usage{
networkDFfromMPoutput(
  chosen.time,
  MPoutput,
  groupNames = NULL,
  sumOverPaths = TRUE,
  sumOverStrains = TRUE,
  convertToMoles = TRUE
)
}
\arguments{
\item{chosen.time}{the time you want to plot}

\item{MPoutput}{the output from microPopModel()}

\item{groupNames}{Default is NULL which plots all the microbes. To plot a subset of all the groups, specify a vector of strings of the names of the groups you want to plot.}

\item{sumOverPaths}{Logical. Default is TRUE which sums flows between the same nodes even if they are on different metabolic paths}

\item{sumOverStrains}{Logical. Default is TRUE which means the strains are put into their functional group nodes and the flow are summed. When it is FALSE, each strain will have its own node.}

\item{convertToMoles}{Logical. Default is TRUE}
}
\value{
a list containing the edges and nodes
}
\description{
make node and edge data frames from microPop output to use in visNetwork
}
