##' @title Plot p-values of allele count differences between populations.
##' 
##' @description Plots p-values of allele count differences between populations.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param data an object of the \code{\link[GenABEL]{gwaa.data-class}}, 
##' @param allele.cnt an object of the \code{\link[pac.result]{pac.result-class}} class,
##' @param plot.LD a logical indicating whether LD pattern should be plotted on top of p-values (experimental),
##' @param legend.pos position of a legend. If null, no legend is shown. If 'default', the default settings are used,
##' @param ...
##' @details Plots p-values of allele count differences between populations.
##' @return null 
##' @examples
##'  \dontrun{
##'  }
##' @keywords plot, allele counts, populations
##' @export
plot.pac <- function(data, allele.cnt, plot.LD=F, legend.pos='default', ...) { 
  chromosomes <- sort(unique(as.numeric(levels(data@gtdata@chromosome))))
  if (length(chromosomes) > 1) {
    if (plot.LD) {
      warning("More than one chromosome supplied. Cannot plot LD pattern.")
      plot.LD <- F
    }
    chr.mid <- c()
    for (chr in chromosomes) {
      coord <- which(as.numeric(data@gtdata@chromosome) == chr)
      chr.mid <- get.chr.midpoints(data)
    }
    color <- (as.numeric(data@gtdata@chromosome)%% 2) + 1
    color[color == 1] <- "slateblue"
    color[color == 2] <- "grey"
    pvals <- allele.cnt@p.values
    plot(-log10(pvals), type='h', col=color, xlab="Chromosome", ylab=expression(paste(-log[10],"(p-value)")), xaxt='n', bty='n', las=2, ...)
    axis(1, at=chr.mid, labels=unique(chromosomes))
  } else {
    range <- max(data@gtdata@map) - min(data@gtdata@map) 
    if (range < 1e6) {
      divisor <- 1e3
      prefix <- "kbp"
    } else {
      divisor <- 1e6
      prefix <- "Mbp"
    }
    axis.start <- min(data@gtdata@map)
    axis.stop <- max(data@gtdata@map)
    pvals <- allele.cnt@p.values
    index.snp <- which(allele.cnt@p.values == min(allele.cnt@p.values, na.rm = T ))[1]
    # Do plotting
    max.y <- max(-log10(pvals), na.rm = T)
    plot(-log10(pvals), type='h', col="darkgrey", xlab=prefix, 
         ylab=expression(paste(-log[10],"(p-value)")), xaxt='n', 
         bty='n', las=2, ylim = range(pretty(c(0, max.y))), 
         panel.first = grid(), cex.axis=.8, ...)
    # Plot LD-colored points
    if (plot.LD) {
      colors <- get.LD.colors(data, chr=chromosomes[1], index.snp)
      points(-log10(pvals), col=as.character(colors[[1]]), pch=as.numeric(colors[[2]]), cex=.7)    
    }
    tmp <- seq(0,length(data@gtdata@map), by=100)
    axis(1, at=tmp, cex.axis=.8,
         labels=round(seq(floor(axis.start/divisor),ceiling(axis.stop/divisor), along.with=tmp),digits=2))
    
    #plot legend
    if (!is.null(legend.pos)) {
      if (legend.pos == 'default') {
        legend.pos <- c(tmp[length(tmp)]-150, max.y)
      }  
      legend(legend.pos, legend=c("(0.8-1.0]","(0.6-0.8]", "(0.4-0.6]", "(0.2-0.4]", "[0.0-0.2]"), 
             pch=19, bty='n', 
             col=c('lightsalmon','lightpink1','lightblue2',"aquamarine2","bisque2"), 
             cex=.8, title=expression(r^2), y.intersp = 0.8)
    }
  }
}
