##'@title SNP_box_twoWay
##'@description Same as SNP_box but based on two SNPs. 
##'This function makes box plots, showing the Genotype x Phenotype map in 2. 
##'Works for both outbreed (nine boxes) and inbreed (four boxes) data. 
##'Can be used to look for epistasis.
##'
##'@param GWASdata Object of class gwaa.data
##'@param marker1 The name or index of one marker/SNP in GWASdata
##'@param marker2 The name or index of a second marker/SNP in GWASdata
##'@param trait The name of a phenotype in GWASdata
##'@param legend - If true, a legend will be added showing the alternative genotypes of the two markers/SNPs
##'@param legendPos - A string with a keyword, specifying the position of the legend. Default = “topleft"
##'@param returnGeno - If true, the genotype matrix with the two markers will be returned
##'@param returnBox - If true, the output-object from the function boxplot will be returned
##'@return boxes
##'@author Simon Forsberg  <\email{simon.forsberg@slu.se}>
##'


SNP.box.twoWay <- function(data, marker1, marker2, trait, legend = F, legendPos = "topleft",returnGeno = F, returnBox = F, ...){
  marker1.geno <- as.genotype.gwaa.data(data[,marker1])
  marker2.geno <- as.genotype.gwaa.data(data[,marker2])
  
  y <- phdata(data)[,trait]
  box <- boxplot(y ~ as.matrix(marker1.geno)*as.matrix(marker2.geno))
  axis(1, at=1:length(box$n), labels=box$n, line=2, lty=0)
  
  if(legend){
    marker1.uniqGeno <- paste(unique(marker1.geno[,1]), collapse=" ")
    marker2.uniqGeno <- paste(unique(marker2.geno[,1]), collapse=" ")
    legend(legendPos, c(paste(marker1, ":", marker1.uniqGeno, collapse=""), paste(marker2, ":", marker2.uniqGeno, collapse="")), cex=.8)
  }
  
  if(returnGeno){
    return(list(marker1.geno=marker1.geno, marker2.geno=marker2.geno))
  }
  if(returnBox){
    return(box)
  }
}
