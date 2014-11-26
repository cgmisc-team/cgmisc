##'@title Get middle point for each chromosome
##' 
##' @description Returns a vector with chromosome mid-point coordinates
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param data a \code{\link[GenABEL]{gwaa.data-class}} object produced by GenABEL
##' @return a vector of mid-point coordinates
##' @details Function result is independent of coordinates used. Function is primarily used 
##' for plotting chromosome labels on the x-axis.
##' @examples
##'  \dontrun{
##'  midpoints <- get.chr.midpoints(data)
##'  }
##' @keywords chromosome, midpoint
##' @return vector of chromosomes midpoints
##' @export get.chr.midpoints
get.chr.midpoints <- function(data) {
  chromosomes <- unique(as.numeric(data@gtdata@chromosome))
  if (length(chromosomes) > 0) {
    chr.mid <- c()
    for (chr in chromosomes) {
      coord <- which(as.numeric(data@gtdata@chromosome) == chr)
      chr.mid <- c(chr.mid, floor((coord[1] + coord[length(coord)])/2))
    }
  }
  else {
    stop("No valid chromosome supplied!")
  }
  return(chr.mid)
}