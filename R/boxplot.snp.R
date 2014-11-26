##'@title Box plot with genotypes
##'
##'@description Makes boxplots of the individuals in every genotypic class, based on a SNP.
##' Works for both outbreed (three boxes) and inbreed (two boxes) data. 
##'@author Simon Forsberg <\email{simon.forsberg@@slu.se}>
##'@param GWASdata a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}}
##'@param marker SNP name
##'@param trait trait name
##'@param recode if genotypes are to be recoded to 'AA','Aa','aa'.
##'@return null
##'@export 
boxplot.snp  <- function(GWASdata, marker, trait, recode = F, ...){
  if(recode){
    geno <- as.double.gwaa.data(GWASdata[,marker])
    geno[geno[,1] == 0,1] <- "AA"
    geno[geno[,1] == 1,1] <- "Aa"
    geno[geno[,1] == 2,1] <- "aa"
  }else{
    geno <- as.genotype.gwaa.data(GWASdata[,marker])
  }
  y <- phdata(GWASdata)[,trait]
  box <- boxplot(y ~ as.matrix(geno), col = c('lightblue2','lightpink1','lightsalmon'))
  axis(1, at=1:length(box$n), labels=paste("n =", box$n), line=2, lty=0)
  genotypes <- as.character(unique(geno[,1]))
  if(length(genotypes) == 1){
    axis(1, at=1, labels=genotypes)
  }
}
