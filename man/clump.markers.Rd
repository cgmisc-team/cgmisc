% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/clump.markers.r
\name{clump.markers}
\alias{clump.markers}
\title{Clump markers according to their LD.}
\usage{
clump.markers(data, gwas.result, chr = 1, bp.dist = 250000, p1 = 1e-04,
  p2 = 0.01, r2 = 0.5, image = F, verbose = F)
}
\arguments{
\item{data}{data object in \code{\link[GenABEL]{gwaa.data-class}},}

\item{gwas.result}{\code{\link[GenABEL]{gwaa.scan-class}} object with association test results,}

\item{chr}{chromosome to be clumped,}

\item{bp.dist}{threshold for inter-marker distance,}

\item{p1}{threshold for index markers,}

\item{p2}{threshold for clumping,}

\item{r2}{threshold for LD,}

\item{image}{a logical indicating whether to plot clumping results or not,}

\item{verbose}{a logical indicating whether to print clumping results as it proceeds}
}
\value{
a list of clumps
}
\description{
\code{clumpMarkers} implements clumping procedure (as described in PLINK documentation) on
a \code{\link[GenABEL]{gwaa.data-class}} object.
}
\examples{
\dontrun{
   clumps <- clump.markers(data.qc0, gwas.result = an0, chr = 6, bp.dist = 250e3, p1 = 0.0001, p2 = 0.01, r2 = 0.5, image=T)
}
}
\author{
Marcin Kierczak <\email{Marcin.Kierczak@imbim.uu.se}>
}
\references{
\url{http://pngu.mgh.harvard.edu/~purcell/plink/clump.shtml}
}

