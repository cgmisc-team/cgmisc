##'@title Calculates heterozygosity for windows
##'@description Calculates average heterozygosity for overlapping windows 
##'produced by \code{get.overlapping.windows} function
##'@author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>, Jagoda Jablonska <\email{jagoda100jablonska@@gmail.com}>
##'@param data data object in \code{\link[GenABEL]{gwaa.data-class}}
##'@param LW a list of windows and coordinates produced by \code{get.overlapping.windows} function.
##'@param progress logical if the progress bar is to be displayed
##'@return Matrix of average heterozygosities for every window
##'@seealso \code{\link[cgmisc]{get.overlapping.windows}}
##'@keywords heterozygosity, window 
##'
##'@export het.overlap.wind

het.overlap.wind <- function (data, LW, progress = TRUE){
  LW <- data.frame(LW[2])
  heter.zyg <- matrix(nrow = (nrow(LW)), ncol = 2)
  heter.zyg[,1]  <- seq(nrow(LW))
  gtypes <- data@gtdata
  if (progress) {
    pb <- txtProgressBar(min=0, max=nrow(LW), initial=0, style=3)
  }
  for(i in seq(nrow(LW))){
    if (progress) {
     setTxtProgressBar(pb=pb, value=i)
    }
    markers <- which(LW[i,] == TRUE, arr.ind= T)
    markers <- as.vector(markers[,2])
    my.gtps <- summary(gtypes[,markers])
    p  <- 1 - my.gtps$Q.2
    q  <- my.gtps$Q.2
    het.tmp  <- 1 - (p**2 + q**2)
    tmp  <- 1/length(het.tmp)*sum(p**2+q**2)
    het.window  <- 1 - tmp
    heter.zyg[i,2]  <-  het.window
}
return (heter.zyg)
}