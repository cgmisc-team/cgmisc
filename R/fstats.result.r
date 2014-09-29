#' @title The \code{fstats.result} class
#'
#' This class contains results returned by the \code{\link[cgmisc]{compute.Fstats}} function.
#'
#' @section Slots: 
#'  \describe{
#'    \item{\code{pops}:}{A \code{list} of per-population statistics}
#'    \item{\code{glob}:}{A \code{data.frame} of global statistics}
#'  }
#'
#' @name fstats.result
#' @rdname fstats.result
#' @aliases fstats.result-class 
#' @exportClass fstats.result
#' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
setClass("fstats.result", representation(pops = "list", glob = "data.frame"))
