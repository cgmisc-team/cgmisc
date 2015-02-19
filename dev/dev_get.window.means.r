# x - a vector of some values per each marker
# LW - a matrix of logicals returned by get.overlapping.windows
get.window.means <- function(x, LW) {
  X <- matrix(rep(x, times=nrow(LW)), nrow=nrow(LW), byrow=T)
  # The trick is to set all the values not-belonging to a window to NA and tell the mean() function
  # to skip the NAs
  X[LW == F] <- NA
  means <- apply(X, 1, mean, na.rm=T)
  return(means)
}
