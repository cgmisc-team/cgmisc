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
##' plot.manhattan.genes(data=data.mh1.qc1, gwas.result=mm.mh1, chr=36, region=c(16e6,21e6), 
##'                   index.snp="BICF2P12960", bonferroni=F, legend.pos="default")
##' }
##' @seealso \code{\link[GenABEL]{gwaa.data-class}}, \code{\link[GenABEL]{scan.gwaa-class}}
##' @export plot.manhattan.genes
##'
plot.manhattan.genes <- function(data, gwas.result, chr, region, index.snp, p.value=0.05, mafThreshold=.05, bed.path=NULL) {  
  require(grid)
  require(ggplot2)
  #require(wq)
  #source("src/plot.genes.R")

  # layOut is copied from the old version of the wq package. 
  # The layOut function is no longer provided by the wq. 
  layOut <- function (...) {
    require(grid)
    x <- list(...)
    n <- max(sapply(x, function(x) max(x[[2]]))) 
    p <- max(sapply(x, function(x) max(x[[3]])))    
    pushViewport(viewport(layout = grid.layout(n, p)))
    
    for (i in seq_len(length(x))) {
      print(x[[i]][[1]], vp = viewport(layout.pos.row = x[[i]][[2]], 
                                       layout.pos.col = x[[i]][[3]]))
      
    }
  }
  
  # Get the minor allele frequency for each SNP
  getMAF <- function(data, region) {
    summ <- summary(data[,region])
    tmp_maf <- (2 * summ$P.22 + summ$P.12) / (2 * summ$NoMeasured)
    tmp_ma_count <- summ$P.22 + 0.5 * summ$P.12
    tmp_het <- summ$P.12 / summ$NoMeasured
    maf <- data.frame(snp=as.character(row.names(summ)), 
                      chr=summ$Chromosome, pos=data@map[region], 
                      allele=summ$A2, maf=tmp_maf, maCnt=tmp_ma_count, heterozygosity=tmp_het)
  }

  mafThreshold <<- mafThreshold
  
  # Get the markers in the specified region and categorize them according to their p-value
  startCoord <- region[1]
  stopCoord  <- region[2]
  myChromosome <- data@gtdata[ ,which(data@gtdata@chromosome == chr)]
  area <- which(myChromosome@map >= startCoord & myChromosome@map <= stopCoord)
  r2matrix <- r2fast(myChromosome, snpsubset = area)
  r2matrix[lower.tri(r2matrix)] <- t(r2matrix)[lower.tri(r2matrix)]
  markers <- which(data@gtdata@snpnames %in% names(area)) 
  markers.coords <<- data@gtdata@map[markers]
  pvals <- -log10(gwas.result@results$P1df[markers])

  
  r2vec <- r2matrix[index.snp,]
  r2vec[is.na(r2vec)] <- -1
  r2col <<- cut(r2vec, breaks=c(1.0,0.8,0.6,0.4,0.2,0.0,-1), 
               labels=rev(c("(0.8-1.0]","(0.6-0.8]", "(0.4-0.6]", "(0.2-0.4]", "[0.0-0.2]", "INDEX")), include.lowest=T)
  
  df2 <<- data.frame(markers = markers.coords, pvals = pvals)
  maf <- getMAF(myChromosome, area)
  maf2 <<- pmin(maf$maf, (1 - maf$maf))
  gwas_palette <- rev(c("#9E0508","tomato","chartreuse3","cyan3","navy","black"))

  # First plot the manhattan plot
  manhattan_plot <- ggplot(df2, aes(x=markers.coords, y=pvals, colour = r2col)) +
    coord_cartesian(xlim = region) +
    scale_colour_manual(name=expression(r^2),
                        values=gwas_palette,
                        labels = c("INDEX","[0.0-0.2]", "(0.2-0.4]","(0.4-0.6]","(0.6-0.8]","(0.8-1.0]")) +
    scale_shape_manual(name = expression(r^2),
                       values=gwas_palette,
                       labels = c(112, 15, 15, 15, 15, 15)) +
    geom_point(size=2) +
    geom_point(data=df2[index.snp,], aes_string(x=df2[index.snp,1], y=df2[index.snp,2]),
               colour="black", shape=21, fill="white", size=2, guide = 'legend') +
    ylab(expression(-log[10](p-value))) +
    theme_bw() +
    theme(plot.margin=unit(c(5,7,0,5),"mm"), legend.position=c(0.9,0.85),
          legend.key=element_blank(), panel.border=element_blank(),
          legend.key.height = unit(4, "mm")) + 
    theme(axis.line = element_line(color = 'grey'), axis.text.x=element_blank(),
          axis.title.x=element_blank())
    
  
  # Second plot showing the minor allele frequency
  maf_plot <- ggplot(df2, aes(x=markers.coords, y=(maf2))) +
    coord_cartesian(xlim = region) +
    geom_line(color = "#2C7FB8") +
    ylab("MAF") +
    geom_hline(aes(yintercept=mafThreshold), color = "tomato") +
    theme_bw() +
    theme(plot.margin=unit(c(0,7,0,4),"mm"), legend.position="none", panel.border=element_blank(),
          legend.text=element_text(size=2), legend.background = element_blank()) +
    theme(axis.line = element_line(color = 'grey'), axis.text.x=element_blank(),
          axis.title.x=element_blank())
    
  
  # Third plot showing the genes in the specified region
  gene_plot <- plot.genes(region, chr, bed.path=bed.path)
  layOut(list(manhattan_plot, 1:4, 1),
         list(maf_plot, 5, 1),
         list(gene_plot, 6:7, 1))
  
}
