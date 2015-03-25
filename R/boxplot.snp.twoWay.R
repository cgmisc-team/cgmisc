##'@title SNP_box_twoWay
##'@description Same as boxplot.snp but based on two SNPs. 
##'This function makes box plots, showing the Genotype x Phenotype map for two loci. 
##'Works for both outbreed (nine boxes) and inbreed (four boxes) data. 
##'Can be used to look for epistasis.
##'
##'@param data a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}}
##'@param marker1 The name or index of one marker/SNP in data
##'@param marker2 The name or index of a second marker/SNP in data
##'@param trait The name of a phenotype in data
##'@param legend If TRUE, a legend showing the alternative genotypes for the two SNPs will be added to the plot.
##'@return null
##'@author Simon Forsberg  <\email{simon.forsberg@slu.se}>
##'
boxplot.snp.twoWay <- function(data, marker1, marker2, trait, legend = F, ...){
  marker1.geno <- as.genotype.gwaa.data(data[,marker1])
  marker2.geno <- as.genotype.gwaa.data(data[,marker2])
  
  y <- phdata(data)[,trait]
  box <- boxplot(y ~ as.matrix(marker1.geno)*as.matrix(marker2.geno), col = 'lightblue', las = 1, frame=F, cex.axis = .8, xaxt = "n", ...)
  grid()
  axis(1, at=1:length(box$n), labels=paste("n =", box$n), line=2, lty=0)
  
  labels.marker1 <- gsub(pattern = "(.*)\\..*", replacement = "\\1", x = box$names)
  labels.marker2 <- gsub(pattern = ".*\\.(.*)", replacement = "\\1", x = box$names)
  labels.marker2.nrPerClass <- table(labels.marker2)[1]
  tmp <- 1 + (labels.marker2.nrPerClass - 1)/2
  tmp2 <- sapply(X = 2:length(unique(labels.marker2)) - 1, FUN = function(x){tmp + x*labels.marker2.nrPerClass })
  labels.marker2.at <- c(tmp, tmp2)
  axis(1, at = 1:length(labels.marker1), labels = labels.marker1, cex.axis = .8)
  
  labels.marker2.table <- table(labels.marker2)
  j <- labels.marker2.table[1]
  for(i in 2:length(labels.marker2.table)){
    abline(v = j + .5)
    j <- j + labels.marker2.table[1]
  }
  
  if(legend){
    marker1.uniqGeno <- paste(unique(marker1.geno[,1]), collapse=" ")
    marker2.uniqGeno <- paste(unique(marker2.geno[,1]), collapse=" ")
    legend("topleft", c(paste(marker2, ":", marker2.uniqGeno, collapse=""), paste(marker1, ":", marker1.uniqGeno, collapse="")), cex=.8, text.col = c("red", "black"))
    axis(3, at = labels.marker2.at, labels=unique(labels.marker2), cex.axis = .8, col.axis = "red")
  }
  else{
    axis(3, at = labels.marker2.at, labels=unique(labels.marker2), cex.axis = .8)
  }
}