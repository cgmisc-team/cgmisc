##' @title Converts gwaa-data-class GenABEL object to a bigRR-readable format.
##' 
##' @description \code{gwaa2bigRR} converts gwaa-data-class GenABEL object to a bigRR-readable format. 
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@slu.se}> 
##' @param formula the part of formula for fixed effects (the ones with no shrinkage)
##' @param data data object in \code{\link[GenABEL]{gwaa.data-class}},
##' @param trait trait name
##' 
##' @return an object of a bigRR-data.class
##' 
##' @references \url{}
##' 
##' @examples \dontrun{
##'    bigRR.data <- gwaa2bigRR(~sex, data.qc0, gt)
##'} 
##'    
##' @export gwaa.to.bigrr gwaa2bigrr gwaa.to.bigRR gwaa2bigRR
##' @aliases gwaa2bigrr gwaa.to.bigRR gwaa2bigRR

gwaa.to.bigrr <- function(formula, data, trait) {
  # Copyright Marcin.Kierczak@slu.se<mailto:Marcin.Kierczak@slu.se>
  # v. 1.1b, 07.11.2013
  Z1 <- as.matrix(as.numeric(data@gtdata))
  # Imputation
  for (i in 1:dim(Z1)[2]) {
    NAs <- which(is.na(Z1[,i]))
    if (length(NAs) > 0) {
      warning(paste("NAs found in marker", data@gtdata@snpnames[i], ". Imputing values!!!"))
    }
    nonNAs <- which(!is.na(Z1[,i]))
    if (length(nonNAs) == 0) {
      stop(paste("A marker with only NAs found:", data@gtdata@snpnames[i],"!"))
    }
    samp <- sample(nonNAs, length(NAs), replace=T)
    Z1[NAs, i] <- Z1[samp, i]
  }
  Z <- scale(Z1)
  X <- model.matrix(formula, data@phdata)
  y <- data@phdata[,trait]
  # Remove missing genotypes
  missing <- which(is.na(y))
  if (length(missing) > 0) {
    warning(paste("Missing phenotype for individuals: ", data@phdata$id[missing], ". Removing rows!"))
    y <- y[-missing]; X <- X[-missing,]; Z <- Z[-missing,]
  }
  #setClass("bigRR.data", representation(y = "vector", X = "matrix", Z = "matrix"))
  result <- new("bigRR.data", y=y, X=X, Z=Z)
}
gwaa2bigrr <- gwaa.to.bigrr
gwaa.to.bigRR <- gwaa.to.bigrr
gwaa2bigRR <- gwaa.to.bigrr