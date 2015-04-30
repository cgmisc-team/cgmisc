##' @title Plot results of \code{compute.fstats}
##'
##' @description Plot results obtained using \code{\link[cgmisc]{compute.fstats}} function.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param data an object of the \code{\link[GenABEL]{gwaa.data-class}}
##' @param fstats an object of the \code{\link[cgmisc]{fstats.result-class}}
##' @return NULL
##' @examples
##' \dontrun{plot.fstats(data, fst)}
##' @keywords plot FST, fst, Fst, fstats, populations, subpopulations
##' @seealso \code{\link[cgmisc]{compute.fstats}}, \code{\link[cgmisc]{fstats.result-class}}
##' @export plot.Fstats
plot.fstats <- function(data, fstats, ...) {
  chr.mid <- c()
  for (chr in unique(as.numeric(data@gtdata@chromosome))) {
    coord <- which(as.numeric(data@gtdata@chromosome) == chr)
    chr.mid <- get.chr.midpoints(data)
  }
  myfst <- fstats@glob$FST
  color <- 'lightsalmon'

  range <- max(data@gtdata@map) - min(data@gtdata@map) 
  if (range < 1e6) {
    divisor <- 1e3
    prefix <- "kbp"
  } else {
    divisor <- 1e6
    prefix <- "Mbp"
  }

  max.y <- max(myfst, na.rm = T)
  plot(myfst, type='h', col=color, xlab=prefix, ylab=expression(F[ST]),
       xaxt='n', las=1, bty='n', ylim = range(pretty(c(0, max.y))),  
       panel.first=grid(), cex.axis=.8)
  
  axis.start <- min(data@gtdata@map)
  axis.stop <- max(data@gtdata@map)
  
  tmp <- seq(0,length(data@gtdata@map), by=100)
  
  axis(1, at=tmp, cex.axis=.8,
       labels=round(seq(axis.start,axis.stop, along.with=tmp)/divisor,digits=2))
}

plot.Fst <- plot.fstats
