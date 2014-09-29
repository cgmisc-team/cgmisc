#' @title The \code{bigRR.data} class
#'
#' This class contains results of running the \code{\link[cgmisc]{gwaa2bigRR}} function.
#' 
#' @section Slots: 
#'  \describe{
#'    \item{\code{y}:}{a \code{vector} with response variable,}
#'    \item{\code{X}:}{design \code{matrix} for fixed effects,}
#'    \item{\code{Z}:}{design \code{matrix} for random effects.}
#'  }
#'
#' @name bigRR.data
#' @rdname bigRR.data
#' @aliases bigRR.data-class
#' @exportClass bigRR.data
#' @author Marcin Kierczak <\email{Marcin.Kierczak@@slu.se}>
setClass("bigRR.data", representation(y = "vector", X = "matrix", Z = "matrix"))
