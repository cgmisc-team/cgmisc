##' @title Plot LD pattern and MAF in a Manhattan plot
##' 
##' @description Function for plotting local LD pattern (relative to a pre-selected marker) on a 
##' Manhattan plot resulting from genome-wide association study (GWAS). Each marker on the plot is 
##' colored according to its LD with the reference marker. Color codes are not continuous, but LD is 
##' discretized into LD intervals. This is useful for examining signals found in a GWAS study. 
##' In addition, minor allele frequency will be plotted in the lower panel.
##' @param data a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}}
##' @param gwas.result a scan.gwaa class object as used by \code{\link[GenABEL]{scan.gwaa-class}}
##' @param chr chromosome number (name) to be displayed
##' @param region a vector of two coordinates to display
##' @param index.snp index marker of the reference SNP (name of the marker)
##' @param p.value p-value threshold to visualize using red dashed line
##' @param bonferroni logical indicating whether Bonferroni-correction shall be used for the displayed p-value threshold
##' @param mafThreshold  threshold value for minor allele frequency (MAF) -- displayed as a red line on the lower panel
##' @param legend.pos a string specifying position of the legend on the plot: "def.left", "def.right" or any other accepted by \code{\link{legend}}
##' @return NULL
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @keywords Manhattan LD linkage disequilibrium genetics
##' @examples 
##' \dontrun{
##' plot.manhattan.LD(data=data.mh1.qc1, gwas.result=mm.mh1, chr=36, region=c(16e6,21e6), 
##'                   index.snp="BICF2P12960", bonferroni=F, legend.pos="default")
##' }
##' @seealso \code{\link[GenABEL]{gwaa.data-class}}, \code{\link[GenABEL]{scan.gwaa-class}}
##' @export plot.manhattan.ld plot.manhattan.LD
##' @aliases plot.manhattan.LD
##' 
plot.manhattan.ld <- function(data, gwas.result, chr, region, index.snp, p.value = 0.05, bonferroni = T, mafThreshold = .05, legend.pos = 'def.left') {  
  opar <- par()
  par(las = 2)
  shift <- 2.3
  topMargin <- 0
  # Helper function for retrieving MAF in a region (chromosome)
  getMAF <- function(data, region) {
    summ <- summary(data[,region])
    tmp_maf <- (2 * summ$P.22 + summ$P.12) / (2 * summ$NoMeasured)
    tmp_ma_count <- summ$P.22 + 0.5 * summ$P.12
    tmp_het <- summ$P.12 / summ$NoMeasured
    maf <- data.frame(snp=as.character(row.names(summ)), 
                      chr=summ$Chromosome, pos=data@map[region], 
                      allele=summ$A2, maf=tmp_maf, maCnt=tmp_ma_count, heterozygosity=tmp_het)
  }
  startCoord <- region[1]
  stopCoord  <- region[2]
  myChromosome <- data@gtdata[ ,which(data@gtdata@chromosome == chr)]
  region <- which(myChromosome@map >= startCoord & myChromosome@map <= stopCoord)
  r2matrix <- r2fast(myChromosome, snpsubset = region)
  r2matrix[lower.tri(r2matrix)] <- t(r2matrix)[lower.tri(r2matrix)]
  markers <- which(data@gtdata@snpnames %in% names(region)) 
  markers.coords <- data@gtdata@map[markers]
  idx.marker <- which(data@gtdata@snpnames == index.snp)
  idx.marker.coords <- data@gtdata@map[idx.marker]
  pvals <- -log10(gwas.result@results$P1df[markers])
  plot(markers.coords, pvals, type='n', xlab="Position (Mb)", 
       ylab=expression(-log[10](p-value)), 
       ylim=c(-shift, max(pvals) + 3), 
       axes=F, panel.first=grid(), las=2)
  # Plot grid
  #abline(h=seq(.5, max(pvals) + topMargin, .5), col="grey", lty=3)
  #points(idx.marker.coords, -log10(gwas.mm@results$P1df[idx.marker]), col="red", pch=19, cex=1)
  r2vec <- r2matrix[index.snp,]
  r2vec[is.na(r2vec)] <- -1
  r2col <- cut(r2vec, breaks=c(1.0,0.8,0.6,0.4,0.2,0.0,-1), 
               labels=rev(c("#9E0508","tomato","chartreuse3","cyan3","navy", "black")), include.lowest=T)
  r2pch <- rep(19, length(r2col))
  r2pch[which(r2col == "black")] <- 1
  points(markers.coords, -log10(gwas.result@results$P1df[markers]), col=as.character(r2col), pch=r2pch, cex=.8)
  if (bonferroni) {
    p.value <- -log10(p.value/nsnps(data))
    abline(h=p.value, col="red", lty=2)
  }
  
  # Legend
  if (legend.pos[1] == "def.left") {
    legend(startCoord + 10, max(pvals) + 1, legend=c("(0.8-1.0]","(0.6-0.8]", "(0.4-0.6]", "(0.2-0.4]", "[0.0-0.2]"), pch=15, bty='n', 
           col=c("#9E0508","tomato","chartreuse3","cyan3","navy"), cex=.7, title=expression(r^2))
  } else if (legend.pos[1] == "def.right") {
    legend(stopCoord - 30, max(pvals) + 1, legend=c("(0.8-1.0]","(0.6-0.8]", "(0.4-0.6]", "(0.2-0.4]", "[0.0-0.2]"), pch=15, bty='n', 
           col=c("#9E0508","tomato","chartreuse3","cyan3","navy"), cex=.7, title=expression(r^2))    
  } else {
    if (length(legend.pos) > 1) {
        legend(legend.pos[1], legend.pos[2], legend=c("(0.8-1.0]","(0.6-0.8]", "(0.4-0.6]", "(0.2-0.4]", "[0.0-0.2]"), pch=15, bty='n', 
               col=c("#9E0508","tomato","chartreuse3","cyan3","navy"), cex=.7, title=expression(r^2))    
    } else {
        legend(legend.pos, legend=c("(0.8-1.0]","(0.6-0.8]", "(0.4-0.6]", "(0.2-0.4]", "[0.0-0.2]"), pch=15, bty='n', 
               col=c("#9E0508","tomato","chartreuse3","cyan3","navy"), cex=.7, title=expression(r^2))    
    }
  }
  
  # Plot MAF below Manhattan (shift units lower)
  maf <- getMAF(myChromosome, region)
  maf2 <- pmin(maf$maf, (1 - maf$maf))
  lines(markers.coords, (4 * maf2) - shift, col="#2C7FB8")
  # Plot horizontal line at MAF = 0.05
  abline(h=mafThreshold * 4 - shift, col=rgb(1,0,0,1), lty=2)
  # Plot axes
  step <- (stopCoord - startCoord) / 5
  axis(1, at = seq(startCoord, stopCoord, by=step), labels=format(seq(startCoord, stopCoord, by=step)/1e6, scientific=F, digits=3), las=1)
  axis(2, at = 0:(max(pvals) + topMargin + 1), las=2)
  axis(4, at = c(2 - shift, 1 - shift, 0.2 - shift, 0 - shift), labels = c(.5, .25, .05, 0))
  mtext("MAF", side=2, at=0.2-shift, outer=F)
  par(opar)
}
plot.manhattan.LD <- plot.manhattan.ld
