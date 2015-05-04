##'@title Get adjacent markers within a distance around the given marker
##'@description Returns markers adjacent to a given marker (within a preset distance).
##'@author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##'@param data a gwaa.data class object as used by the \code{\link[GenABEL]{gwaa.data-class}},
##'@param marker a name of the central marker,
##'@param size.bp window size in basepairs,
##'@return a matrix of marker adjacency
##'@keywords adjacent markers
##'@examples 
##'\dontrun{
##' adj.markers  <- get.adjacent.markers(data, 'BICF2P647127') 
##' }
##'@export

get.adjacent.markers <- function(data, marker, size.bp=1e3) {
  chr <- data@gtdata@chromosome[marker]
  chr.data <- data[,which(data@gtdata@chromosome == chr)]
  coord <- chr.data@gtdata@map[marker]
  start <- pmax(min(chr.data@gtdata@map), coord - size.bp)
  end <- pmin(coord + size.bp, max(chr.data@gtdata@map))
  tmp <- chr.data[,which(chr.data@gtdata@map >= start & chr.data@gtdata@map <= end)]
  return(as.double(tmp@gtdata))
}
