##' @title Plot marker clumps on Manhattan plot.
##' 
##' @description Plot clumps resulting from running the \code{\link[cgmisc]{clump.markers}} function.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param gwas.result an object of the \code{\link[GenABEL]{gwaa.data-class}}, 
##' @param clumps a result of running the \code{\link[cgmisc]{clump.markers}} function,
##' @param chr chromosome to display,
##' @param region a vector of start and stop coordinates of a region to display,
##' @param clambda a logical indicating whether corrected Pc1df p-values are
##' to be used.
##' @return NULL
##' @examples
##'  \dontrun{plot.clumps(data, myclumps, 1, c(14172, 19239))}
##' @keywords plot clumps clumping
##' @seealso \code{\link[cgmisc]{clump.markers}}
##' @export plot.clumps
plot.clumps <- function(gwas.result, clumps, chr, region, clambda = F) {
  if (length(clumps) > 0) {
    par(mfrow=c(1,1))
    minCoord <- region[1]
    maxCoord <- region[2]
    region <- which(gwas.result@annotation$Chromosome == chr & gwas.result@annotation$Position > minCoord & gwas.result@annotation$Position < maxCoord)
    if (length(region) == 0) {
      stop("Supplied region appears to be a marker desert. Emptiness...")
    }
    if (clambda) {
      pvals <- -log10(gwas.result@results$Pc1df[region])
    } else {
      pvals <- -log10(gwas.result@results$Pc1df[region])
    }
    coords <- gwas.result@annotation$Position[region]
    plot(coords, pvals, pch=19, cex=1, ann = FALSE, xaxt='n', yaxt='n', 
         bty='n', type = 'n', 
         ylim = c(-length(clumps) - 1, max(pvals)), las = 2)
    grid()
    step <- (max(coords) - min(coords)) / 5
    axis(1, at = seq(min(coords), max(coords), by=step),
         labels = format(seq(min(coords), max(coords), by=step)/1e6, 
                       scientific=F, digits=3))
    axis(2, at = seq(0, max(pvals), 2), 
         labels = T, las = 1, cex.axis = 0.6)
    #axis(4, at=seq(-1, -1-length(clumps), -1),labels=abs(seq(-1, -1-length(clumps), -1)))
    mtext('Position (Mb)', 1, 3)
    #mtext(expression(clumps~'            '~-log[10](p-value)), 2, 3)
    mtext(expression(-log[10](p-value)), 2, 3)
    points(coords, pvals, pch=19, cex=.7)
    mycols <- colorRampPalette(colors=c("slateblue","grey"))
    cols <- mycols(length(clumps))
    for (i in 1:length(clumps)) {
      x <- rep(min(clumps[[i]]$coord), times=2)
      y <- c(-1-i, max(pvals))
      #lines(x,y, col="tomato", lty=2)
      lines(c(min(clumps[[i]]$coord), max(clumps[[i]]$coord)), c(-1-i, -1-i), col=cols[i])
      #points(clumps[[i]]$coord, -log10(clumps[[i]]$pval), col=cols[i], type='l', lwd=2)
      points(clumps[[i]]$coord, -log10(clumps[[i]]$pval), col=cols[i], cex=.9, pch=19)
      #abline(v=(min(clumps[[i]]$coord)+max(clumps[[i]]$coord))/2, col="tomato", lty=2)
      #points(clumps[[i]]$coord, rep(-1-i, times=length(clumps[[i]]$coord)), col=cols[i], pch=19, cex=.9)
    }
  }
  else { stop("The list of clumps is empty!") }
}
