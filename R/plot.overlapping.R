##'@title Plots average heterozygosities for overlapping windows
##'
##'@description Plots average heterozygosities for overlapping windows and the distribution of heterozygosity.
##'@author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>, Jagoda Jablonska <\email{jagoda100jablonska@@gmail.com}>
##'@param LW a list of windows and coordinates returned by \code{get.overlapping.windows} function
##'@param heter.zyg a matrix of heterozygosities returned by \code{\link[cgmisc]{het.for.overlap.wind}}
##'@return null
##'@seealso \code{\link[cgmisc]{get.overlapping.windows}}, \code{\link[cgmisc]{het.for.overlap.wind}}
##'@keywords heterozygosity, window 
##'
##'@export

plot.overlapping  <- function(LW, heter.zyg){
  mids  <- as.vector(LW[[1]][[2]])
  het  <- as.vector(heter.zyg[,2])   
  plot(x = mids, y = het, xlab = 'windows mids', ylab ='average heterozygosities',
       type = 'h', main = 'Heterozygosities per windows', col = 'dark red', las=1)
  hist(het, breaks = 50, main = 'Heterozygosities distribution', xlab = 'Heterozygosity',
       col = 'darkolivegreen4', las = 1)
}