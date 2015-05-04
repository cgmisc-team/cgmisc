##' @title Plot endogenous retroviral sequences in canine genome.
##' 
##' @description Plot endogenous retroviral sequences (ERV) in canine genome.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>, Patric Jern <\email{Patric.Jern@@imbim.uu.se}>
##' @param chr chromosome, e.g. "chr17", 
##' @param coords a vector of coordinates,
##' @param src source database, currently only canFam3 is supported,
##' @param scale scaling factors for x and y axis.
##' @details Plots -log10(Retrotector score) on the y axis (the higher the score is, the more probable that the hit represents an actual ERV)
##' and coordinates in Mbp on the x axis. Width of the bars covers entire ERV. We display >300. 
##' @return null 
##' @examples
##'  \dontrun{plot.erv("chr16", coords=c(1e6,2e6))
##'  }
##' @keywords plot, ERV
##' @export
plot.erv <- function(chr=NA, coords=c(NA,NA), src="canFam3cgmisc.db", scale=c(.2,.1)) {
  result <- get.erv(chr=chr, coords=coords, src=src)
  tmp <- log10(result$score)
  my.ylim <- c(0, max(tmp) + scale[2] * max(tmp))
  my.xlim <- c(min(result$start) - scale[1]*min(result$start), max(result$end) + scale[1]*max(result$end))
  plot(0, xlim = my.xlim, ylim=my.ylim, type='n', bty='n', las=1, xaxs='i', cex.axis=.7, cex.lab=.7, 
       xaxt='n', xlab="Mbp", ylab=expression(-log[10](score)))
  grid()
  p <- pretty(my.xlim)
  axis(1, at = p, labels=p/1e6, cex.axis=.7)
  for (i in 1:dim(result)[1]) {
    r <- result[i,]
    scol <- r$strand
    if (scol == "P") {
      scol <- "slateblue"
    } else {
      scol <- "olivedrab"
    }
    polygon(c(r$start, r$end, r$end, r$start, r$start), 
            c(0,0,log10(r$score),log10(r$score), 0), col = scol, border = 1)
    points(mean(c(r$start, r$end)), log10(r$score), col = scol, pch=15, cex=.8)
  }
  opar <- par()
  par(xpd=TRUE) 
  legend(my.xlim[1],-1.5, col=c("slateblue","olivedrab"), legend=c("+ strand","- strand"), bty = 'n', pch=15, cex=.8) 
  par(opar) 
}
