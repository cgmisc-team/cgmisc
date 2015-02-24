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
boxplot.snp  <- function(data, marker, trait, recode = F, ...){
  if(recode){
    geno <- as.double.gwaa.data(data[,marker])
    geno[geno[,1] == 0,1] <- "AA"
    geno[geno[,1] == 1,1] <- "Aa"
    geno[geno[,1] == 2,1] <- "aa"
  }else{
    geno <- as.genotype.gwaa.data(data[,marker])
  }
  y <- phdata(data)[,trait]
  max.y <- max(y, na.rm=T)
  box <- boxplot(y ~ as.matrix(geno), col = 'lightblue', las = 1, frame=F, cex.axis = .8,  ylim = range(pretty(c(0, max.y))))
  grid()
  boxplot(y ~ as.matrix(geno), col = 'lightblue', las = 1, frame=F, cex.axis = .8,  ylim = range(pretty(c(0, max.y))), add = T)
  axis(1, at=1:length(box$n), labels=paste("n =", box$n), line=2, lty=0)
  axis(1, labels=FALSE)
  genotypes <- as.character(unique(geno[,1]))
  if(length(genotypes) == 1){
    axis(1, at=1, labels=genotypes)
  }
}
