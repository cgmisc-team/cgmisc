##'@title Plots average heterozygosities for overlapping windows
##'
##'@description Plots average heterozygosities for overlapping windows and the distribution of heterozygosity.
##'@author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>, Jagoda Jablonska <\email{jagoda100jablonska@@gmail.com}>
##'@param LW a list of windows and coordinates returned by \code{get.overlap.windows} function
##'@param heter.zyg a matrix of heterozygosities returned by \code{\link[cgmisc]{het.overlap.wind}}
##'@return null
##'@seealso \code{\link[cgmisc]{get.overlap.windows}}, \code{\link[cgmisc]{het.overlap.wind}}
##'@keywords heterozygosity, window 
##'
##'@export plot.overlap

plot.overlap  <- function(LW, heter.zyg, ...){
  mids  <- as.vector(LW[[1]][[2]])
  het  <- as.vector(heter.zyg[,2])  
  max.y <- max(het, na.rm = T)
  plot(x = 1:length(het), y = het, xlab = 'window', ylab ='average heterozygosities',
       type = 'h',lwd=2, col = 'darkseagreen3', las=1, frame = F, ylim = range(pretty(c(0, max.y))), axes = F, ...)
  axis(1, at = seq(1,length(het), length.out = 100), labels = F)
  axis(1,lwd = 0)
  axis(2, las = 1)

}