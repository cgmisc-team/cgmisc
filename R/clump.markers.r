##' @title Clump markers according to their LD.
##' @description \code{clumpMarkers} implements clumping procedure (as described in PLINK documentation) on
##' a \code{\link[GenABEL]{gwaa.data-class}} object.  
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}> 
##' @param p1 threshold for index markers,
##' @param p2 threshold for clumping,
##' @param r2 threshold for LD,
##' @param bp.dist threshold for inter-marker distance,
##' @param chr chromosome to be clumped,
##' @param gwas.result \code{\link[GenABEL]{gwaa.scan-class}} object with association test results,
##' @param data data object in \code{\link[GenABEL]{gwaa.data-class}},
##' @param image a logical indicating whether to plot clumping results or not,
##' @param verbose a logical indicating whether to print clumping results as it proceeds
##' 
##' @return a list of clumps
##' 
##' @references \url{http://pngu.mgh.harvard.edu/~purcell/plink/clump.shtml}
##' 
##' @examples \dontrun{
##'    clumps <- clump.markers(data.qc0, gwas.result = an0, chr = 6, bp.dist = 250e3, p1 = 0.0001, p2 = 0.01, r2 = 0.5, image=T)
##' } 
##'    
##' @export  
clump.markers <- function(data, gwas.result, chr=1, bp.dist=250e3, p1=0.0001, p2=0.01, r2=0.5, image=F, verbose=F) {
  an  <- gwas.result
  data.chr <- data[,data@gtdata@chromosome == chr]
  result <- gwas.result[gwas.result@annotation$Chromosome == chr,]
  result.sorted <- result[order(result$P1df),]
  signif.p1 <- rownames(result.sorted[result.sorted$P1df <= p1,])
  signif.p2 <- rownames(result[result$P1df <= p2,])
  data.signif <- data.chr[,data.chr@gtdata@snpnames %in% signif.p2]
  r2matrix <- r2fast(data.signif)
  r2matrix[lower.tri(r2matrix)] <- t(r2matrix)[lower.tri(r2matrix)] 
  #image(r2matrix)
  d <- as.matrix(dist(cbind(data.signif@gtdata@map, rep(0, times=length(signif.p2)))))
  #image(d)
  clumpmatrix <- matrix(rep(0, times=length(signif.p2)^2), nrow=length(signif.p2), ncol=length(signif.p2))
  clumpmatrix[which(d <= bp.dist)] <- clumpmatrix[which(d <= bp.dist)] + 1
  clumpmatrix[which(r2matrix >= r2)] <- clumpmatrix[which(r2matrix >= r2)] + 3
  colnames(clumpmatrix) <- colnames(d)
  rownames(clumpmatrix) <- rownames(d)
  #image(clumpmatrix)
  marker.names <- data.signif@gtdata@snpnames
  used <- rep(0, times=length(signif.p2))
  clumps <- list()
  for (i in 1:length(signif.p1)) {
    marker <- signif.p1[i]
    marker.index <- which(signif.p2 == marker)
    if (used[marker.index] == 0) {
      used[marker.index] <- 1
      newClump <- list()
      clump <- which(clumpmatrix[marker,] == 4)
      unused <- which(used[clump] == 0)
      clump <- clump[unused]
      if (length(clump) > 0) {
        p <- paste("Marker ", marker, " clumps with markers: ", paste(marker.names[clump], collapse=', '), sep="")
        snpnames <- c(marker, marker.names[clump])
        newClump[["snpnames"]] <- snpnames
        newClump[["chr"]] <- as.character(an@annotation$Chromosome[which(rownames(an@results) %in% snpnames)])
        newClump[["coord"]] <- an@annotation$Position[which(rownames(an@results) %in% snpnames)]
        newClump[["pval"]] <- an@results$P1df[which(rownames(an@results) %in% snpnames)]
        clumps[[marker]] <- newClump
        if (verbose) {
          print(p)
        }
      }
      #else {warning("No clumps found in dataset!")}
      used[clump] <- 1
    }
  }
  if (image == T) {
    par(mfrow=c(1,3))
    image(r2matrix, col=rev(heat.colors(100)), main="r2 matrix")
    image(d, col=rev(heat.colors(100)), main="distance matrix")
    image(clumpmatrix, col=c("cornsilk1","blue","tomato","red"), main="clumping matrix")
  }
  if (length(clumps) == 0) {
    warning("No clumps found in dataset!")
  }
  
  clumps
}
