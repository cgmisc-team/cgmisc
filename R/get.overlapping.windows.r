##'@title Get overlapping windows
##'@description Divides chromosome into a number of overlapping windows.
##'Windows and overlap size is defined in terms of a fixed number of base pairs. 
##'One may use returned windows for further analyses.
##'@author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##'@param data data object in \code{\link[GenABEL]{gwaa.data-class}},
##'@param chr chromosome to be divided into windows - if not specified, the entire genome is taken into consideration,
##'@param size size of the window,
##'@param overlap size of the overlap.
##'
##'@return List containing windows coordinates and LW matrix with windows and markers.
##'@seealso \code{\link[cgmisc]{het.overlap.wind}}
##'@keywords overlapping, windows 
##'
##'@export get.overlap.windows 

get.overlapping.windows <- function(data, chr, size = 125e4, overlap = 25e2) {
  if(!missing(chr)){
    data<-data[,data@gtdata@chromosome == chr]
  } else {
    data <- data
  }
  if (size > (rev(data@gtdata@map)[1] - data@gtdata@map[1])){
    stop("Size parameter cannot be bigger than size of map!")
  }
  if(overlap > (rev(data@gtdata@map)[1] - data@gtdata@map[1])){
    stop("Size of overlap cannot be bigger than size of map!")
  }
  half.window <- floor(0.5 * size) 
  map <- data@gtdata@map
  centers <- seq(map[1] + half.window, rev(map)[1] - half.window, by=(size - overlap))
  starts <- centers - half.window
  stops <- centers + half.window
  windows <- data.frame(start=starts, midpoint=centers, stop=stops)
  LW <- matrix(FALSE, nrow = length(starts), ncol=length(map))
  for (i in 1:length(starts)) {
    LW[i, (map >= starts[i] & map <= stops[i])] <- TRUE
  }
  return(list(windows, LW))
}

