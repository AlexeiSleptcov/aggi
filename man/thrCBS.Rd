% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thrCBS.R
\name{thrCBS}
\alias{thrCBS}
\title{Detection CNA}
\usage{
thrCBS(data, method = c("MAD", "SD"), thr = NULL, n = 1)
}
\arguments{
\item{data}{CBS data}

\item{method}{method of detection threshold. Default \code{MAD}}

\item{n}{numeric, multiple to method}

\item{the}{numeric, if need to establishe custom threshold}
}
\description{
Detection CNA using MAD or SD.
}

