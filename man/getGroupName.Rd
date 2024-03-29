% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGroupName.R
\name{getGroupName}
\alias{getGroupName}
\title{Convert strain name to its group name
e.g. 'Bacteroides.1' becomes 'Bacteroides'
updated (Dec 2019) so that MFG names can contain dots}
\usage{
getGroupName(xname, microbeNames)
}
\arguments{
\item{xname}{a string (may be strain name or something else)}

\item{microbeNames}{vector of strings of microbial group names}
}
\value{
group name (string) if xname is a strain name. If xname is not a the name of a strain it will simply return xname unchanged.
}
\description{
Convert strain name to its group name
e.g. 'Bacteroides.1' becomes 'Bacteroides'
updated (Dec 2019) so that MFG names can contain dots
}
