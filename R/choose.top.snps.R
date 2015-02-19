##'@title Select the highest r2 SNPs to a given top marker.
##'
##'@description Calculates r2 and selects SNPs with highest r2 to a given top marker.
##'
##'@author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}> Jagoda Jablonska <\email{jagoda100jablonska@@gmail.com}
##'@param data data object in \code{\link[GenABEL]{gwaa.data-class}}
##'@param chr chromosome number (name) to be displayed
##'@param region a vector of two coordinates
##'@param index.snp the reference SNP (name of the marker)
##'
##'@return matrix with highest r2 SNPs with their r2 values and coordinates 
##'@export choose.top.snps
choose.top.snps <- function(data, chr, region, index.snp) {  
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
  r2vec <- r2matrix[index.snp,]
  r2vec[is.na(r2vec)] <- -1
  r2col <- cut(r2vec, breaks=c(1.0,0.8,0.6,0.4,0.2,-0.1,-1), 
               labels=rev(c("a","b","c","d","e","f")), include.lowest=T)
 'INDEX SNP' -> r2matrix[is.na(r2matrix)]
  ret.mark  <- rev(sort(r2matrix[index.snp, which(r2col == 'a' | r2col == 'b' | r2col == 'f' )]))
  snps  <- as.vector(names(ret.mark))
  coords  <- data[,snps]@gtdata@map
  all  <- data.frame(marker=ret.mark, coord=coords) 
  return(all)
}
