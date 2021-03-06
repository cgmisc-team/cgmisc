% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/get.roh.R
\name{get.roh}
\alias{get.roh}
\title{Runs of homozygosity}
\usage{
get.roh(data, chr, LW, hetero.zyg, threshold = 0.23, strict = TRUE)
}
\arguments{
\item{data}{a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}}}

\item{chr}{number of chromosome}

\item{LW}{a matrix of windows coordinares returned by \code{get.overlap.windows} function}

\item{hetero.zyg}{a matrix of average heterozygosities for windows returned by \code{het.overlap.wind} function}

\item{threshold}{threshold of homozygosity. All windows with homozygosity below this limit will be treated as homozygous.}

\item{strict}{if one window with lower homozygosity should be tolerated}
}
\value{
matrix with columns respectively : 1) number of first window of run, 2) coordinate of first window in a run,
3) coordinate of last window in a run, 4) length of run (number of windows)
}
\description{
Identifies runs of homozygosity in overlapping windows across chromosome.
}
\author{
Marcin Kierczak <\email{Marcin.Kierczak@imbim.uu.se}>, Jagoda Jablonska <\email{jagoda100jablonska@gmail.com}>
}
\seealso{
\code{\link[cgmisc]{get.overlap.windows}} \code{\link[cgmisc]{het.overlap.wind}}
}
\keyword{heterozygosity,}
\keyword{homozygosity,}
\keyword{overlapping}
\keyword{window}

