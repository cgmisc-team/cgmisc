##'@title Plot genes provided in a bed file
##'
##'@description Using information provided in a bed file, plot genes, together with their names 
##'and strand information. 
##'@author Veronika Scholz <\email{veronikascholz@@gmail.com}>
##'@param region \code{\link[GenomicRanges]{IRanges}}
##'@param bed.path is a path to a bed file containing gene information.
##'@details When using UCSC-provided bigBed annotation files you need to first use the 
##'bigBedToBed tool (available from: http://hgdownload.cse.ucsc.edu/admin/exe/YOUR_OS_DIRECTORY) to 
##'get the right format. At this point you can also extract the desired chromosome.
##'Then you need to strip the last column of the bed file using, for instance awk: 
##'awk '{out=""; for(i=1;i<=15;i++){out=out" "$i}; print out}' my.bed > myFinal.bed
##'Here, we provide a Broad Improved Canine Annotation v.1 protein coding genes for the canine chromosome 2 
##'(see: https://www.broadinstitute.org/ftp/pub/vgb/dog/trackHub/hub.txt).
##'@return ggplot2 plot
##'@export 
plot.genes <- function(region, chr, bed.path=NULL) {
  targetRanges <- IRanges(region[1], region[2])
  reads <- import.bed(con=bed.path, asRangedData=F)
  gr <- GRanges(seqnames=Rle(paste("chr", chr, sep="")), targetRanges)
  hits <- findOverlaps(reads, gr)
  genes <- reads[queryHits(hits)]
  df <- data.frame(start = start(ranges(genes)), end = end(ranges(genes)),
                   name = elementMetadata(genes)$name, strand = strand(genes))
  df_unique <<- df[!duplicated(df[,3]),]
  region_length <- region[2]-region[1]
  df_unique$ymin <<- 0
  for (i in 1:(length(df_unique$start)-1)) {
    for (j in (i+1):length(df_unique$start)) {
      if ((df_unique$end[i]+(region_length/100)) >= df_unique$start[j] ) {
        df_unique$ymin[j] <<- df_unique$ymin[i]+2
      }
    }
  }
  df_unique$seg_start <<- ifelse(df_unique$strand == "+", df_unique$start, df_unique$end)
  df_unique$seg_end <<- ifelse(df_unique$strand == "+", df_unique$end, df_unique$start)
  
  
  gene_plot <- ggplot(df_unique, aes(xmin = df_unique$start, xmax = df_unique$end,
                                     ymin = df_unique$ymin, ymax = df_unique$ymin + 0.7)) +
    coord_cartesian(xlim = region) +
    scale_x_continuous(breaks = seq(region[1], region[2], 1000000),
                       labels = round(seq(region[1], region[2], 1000000), 2)) +
    
    scale_y_continuous(expand=c(0.1,0.3)) +
    xlab("Position (bp)") +
    geom_rect(alpha =0.2, color = rgb(0,114,178, maxColorValue=256), 
              fill = rgb(0,114,178, maxColorValue=256)) +
    geom_segment(aes(x=df_unique$seg_start, xend= df_unique$seg_end,
                     y= df_unique$ymin+1.0, yend =df_unique$ymin+1.0),
                 arrow=arrow(length = unit(3,"mm")), size=1, alpha=0.7, color = rgb(0,114,178, maxColorValue=256)) +
    geom_text(aes(label = df_unique$name, x = df_unique$start-10000,
                  y = df_unique$ymin + 1.3),
                  size = 2, angle=25) +
    theme_bw() +
    theme(plot.margin=unit(c(0,7,0,12),"mm"), panel.border=element_blank()) +
    theme(axis.line = element_line(color = 'grey'),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank())
  return(gene_plot)
}
