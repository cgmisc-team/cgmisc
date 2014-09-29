##' @title Plot results of \code{compute.fstats}
##' 
##' @description Plot results obtained using \code{\link[cgmisc]{compute.fstats}} function.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param data an object of the \code{\link[GenABEL]{gwaa.data-class}} 
##' @param fstats an object of the \code{\link[cgmisc]{fstats.result-class}} 
##' @return NULL
##' @examples
##'  \dontrun{plot.fstats(data, fst)}
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
  col.fun <- wes.palette(n = 2, name = "GrandBudapest", type = "continuous")
  color <- (as.numeric(data@gtdata@chromosome)%% 2) + 1
  color[color == 1] <- col.fun(2)[1]
  color[color == 2] <- col.fun(2)[2]
  plot(myfst, type='h', col=color, xlab="Chromosome", ylab=expression(F[ST]), 
       xaxt='n', las=1, bty='n', ...)
  axis(1, at=chr.mid, labels=unique(dat@gtdata@chromosome))
}