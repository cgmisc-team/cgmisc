#' @title The \code{pac.result} class
#'
#' This class contains results of running the \code{\link[cgmisc]{pop.allele.counts}} function.
#'
#' @section Slots: 
#'  \describe{
#'    \item{\code{chr}:}{A \code{vector} of chromosome numbers.}
#'    \item{\code{map}:}{A \code{vector} of genomic positions.}
#'    \item{\code{raf.pop1}:}{Reference allele frequency in population 1.}
#'    \item{\code{raf.pop2}:}{Reference allele frequency in in population 2.}
#'    \item{\code{n.pop1}:}{Number of individuals in population 1.}
#'    \item{\code{n.pop2}:}{Number of individuals in population 2.}
#'    \item{\code{p.values}:}{A \code{vector} of p-values.}
#'  }
#'
#' @name pac.result
#' @rdname pac.result
#' @aliases pac.result-class 
#' @exportClass pac.result
#' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
setClass("pac.result", representation(chr = "vector",
                                      map = "vector",
                                      p.values = "vector", 
                                      raf.pop1 = "vector",
                                      raf.pop2 = "vector",
                                      n.pop1 = "vector",
                                      n.pop2 = "vector"))
