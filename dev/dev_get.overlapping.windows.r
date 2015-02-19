dat <- data[,data@gtdata@chromosome == 1]
tmp <- get.overlapping.windows(dat)

get.overlapping.windows <- function(data, size = 125e3, overlap = 25e3) {
  half.window <- floor(0.5 * size) 
  # Get vector of marker coordinates
  map <- data@gtdata@map
  # Retrieve window centers (rev(map)[1] is a trick to get the last element)
  centers <- seq(map[1] + half.window, rev(map)[1] - half.window, by=(size - overlap))
  starts <- centers - window
  stops <- centers + window
  # windows stores the result together with M. M is a P x Q matrix where P is the number of 
  # windows and Q is the number of markers. Initially the LW matrix is filled with FALSE
  # for each window, markers that are inside this window are set to TRUE
  windows <- data.frame(start=starts, midpoint=centers, stop=stops)
  LW <- matrix(FALSE, nrow = length(starts), ncol=length(map))
  for (i in 1:length(starts)) {
    LW[i, (map >= starts[i] & map <= stops[i])] <- TRUE
  }
  return(list(windows, LW))
}
