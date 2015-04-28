plot.manhattan <- function(data, gwas.result, chr, region, index.snp, p.value=0.05, mafThreshold=.05) {  

  library(grid)
  require(ggplot2)
  require(wq)
  source("src/plot.genes.R")

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
  
  startCoord <- region[1]
  stopCoord  <- region[2]
  myChromosome <- data@gtdata[ ,which(data@gtdata@chromosome == chr)]
  area <- which(myChromosome@map >= startCoord & myChromosome@map <= stopCoord)
  r2matrix <- r2fast(myChromosome, snpsubset = area)
  r2matrix[lower.tri(r2matrix)] <- t(r2matrix)[lower.tri(r2matrix)]
  markers <- which(data@gtdata@snpnames %in% names(area)) 
  markers.coords <<- data@gtdata@map[markers]
  idx.marker <- which(data@gtdata@snpnames == index.snp)
  idx.marker.coords <- data@gtdata@map[idx.marker]
  pvals <- -log10(gwas.result@results$P1df[markers])
  
  r2vec <- r2matrix[index.snp,]
  r2vec[is.na(r2vec)] <- -1
  r2col <<- cut(r2vec, breaks=c(1.0,0.8,0.6,0.4,0.2,0.0,-1), 
               labels=rev(c("(0.8-1.0]","(0.6-0.8]", "(0.4-0.6]", "(0.2-0.4]", "[0.0-0.2]", "INDEX")), include.lowest=T)
  
  df2 <- data.frame(markers = markers.coords, pvals = pvals)
  maf <- getMAF(myChromosome, area)
  maf2 <<- pmin(maf$maf, (1 - maf$maf))
  gwas_palette <- rev(c("#9E0508","tomato", "chartreuse3", "cyan3","navy","black" ))
  manhattan_plot <- ggplot(df2, aes(x=markers.coords, y=pvals, colour = r2col)) +
    coord_cartesian(xlim = region) +
    scale_colour_manual(values=gwas_palette, name=expression(r^2)) +
    geom_point(size=3) +
    ylab(expression(-log[10](p-value))) +
    theme_bw() +
    theme(plot.margin=unit(c(5,7,0,5),"mm"), legend.position=c(0.9,0.8), panel.border=element_blank()) +
    theme(axis.line = element_line(color = 'grey'), axis.text.x=element_blank(),
          axis.title.x=element_blank())
  
  maf_plot <- ggplot(df2, aes(x=markers.coords, y=(maf2))) +
    coord_cartesian(xlim = region) +
    geom_line(color = "#2C7FB8") +
    ylab("MAF") +
    geom_hline(aes(yintercept=mafThreshold), color = "tomato") +
    theme_bw() +
    theme(plot.margin=unit(c(0,7,0,4),"mm"), legend.position="none", panel.border=element_blank()) +
    theme(axis.line = element_line(color = 'grey'), axis.text.x=element_blank(),
          axis.title.x=element_blank())
  
  gene_plot <- plot.genes(region, "data/canis_familiaris.prot.bed")
  
  layOut(list(manhattan_plot, 1:4, 1),
         list(maf_plot, 5, 1),
         list(gene_plot, 6:7, 1))
}