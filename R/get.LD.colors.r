##' @title Plot results of \code{get.LD.colors}
##' 
##' @description Returns colour scale based on LD to a chosen index marker.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param data an object of the \code{\link[GenABEL]{gwaa.data-class}} 
##' @param chr number of chromosome of interest
##' @param index.snp numeric index of the index marker 
##' @param region a region for which to obtain LD colors
##' @param wes.par parameters for package wesanderson
##' @return a list with four vectors: 1) LD-based colors, 2) pch values, 3) color scale and 4) breakpoints. 
##' @examples
##'  \dontrun{get.LD.colors(data, chr = 38, index.snp=270) 
##'  }
##' @keywords plot, LD, colours
##' @seealso \code{\link[cgmisc]{plot.pac}}
##' @export
get.LD.colors <- function(data, chr, index.snp, region=NULL, wes.par=c(5, 'Chevalier', 'continuous')) {
  require(wesanderson)
  data  <- data[,data@gtdata@chromosome==chr]
  if (!is.null(region)) {
    data <- data[,which(data@gtdata@map >= region[1] & data@gtdata@map <= region[2])]
  }
  result <- list()
  myBreaks <- c(1.0,0.8,0.6,0.4,0.2,0.0,-1)
  wes.f <- wes.palette(wes.par[1], wes.par[2], wes.par[3])
  myColors <- rev(c(wes.f(5), 'black'))
  r2matrix <- r2fast(data)
  r2matrix[lower.tri(r2matrix)] <- t(r2matrix)[lower.tri(r2matrix)]
  r2vec <- r2matrix[index.snp,]
  r2vec[is.na(r2vec)] <- -1
  r2col <- cut(r2vec, breaks=myBreaks, 
               labels=myColors, include.lowest=T)
  r2pch <- rep(19, length(r2col))
  r2pch[which(r2col == "black")] <- 1
  result[[1]] <-  as.character(r2col)
  result[[2]] <- as.character(r2pch)
  result[[3]] <- myBreaks
  result[[4]] <- as.character(myColors)
  result
}

# Testing function
# myChromosome <- data[,data@gtdata@chromosome==38]
# myChromosome <- myChromosome[,1:501]
# tmp <- get.LD.colors(myChromosome, 251)
# plot(sample(x=c(1:10), size=501, replace=T), pch=tmp[[2]], col=tmp[[1]], cex=.5)
