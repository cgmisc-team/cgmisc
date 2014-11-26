##'@title Compare allele counts in two populations.
##'@description Compare allele counts in two populations using Fisher exact test.
##'@param data a \code{\link[GenABEL]{gwaa.data-class}} object
##'@param pops a vector defining populations. The vector should be as long as 
##'the number of individuals in data and populations are coded as 1 and 2.
##'@return an object of \code{\link[cgmisc]{pac.result-class}}
##'@author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##'@export pop.allele.counts
pop.allele.counts <- function(data, pops, progress=TRUE) {
  
  all <- summary(data@gtdata[,])
  pop.1 <- summary(data@gtdata[which(pops==1),])
  pop.2 <- summary(data@gtdata[which(pops!=1),])
  # Observed allele counts in each population
  obs.1 <- (2 * pop.1$P.11) + pop.1$P.12
  obs.2 <- (2 * pop.2$P.11) + pop.2$P.12 
  
  # Expected allele counts in each population
  freq <- (2 * all$P.11 + all$P.12) / all$NoMeasured # Allele 1 freq in pops jointly
  exp.1 <- 2 * pop.1$NoMeasured * freq
  exp.2 <- 2 * pop.2$NoMeasured * freq
  
  N <- length(freq)
  pvals <- as.numeric(rep(NA, times=N))
  if (progress) {
    pb <- txtProgressBar(min=0, max=N, initial=0, style=3)
  }
  for (i in 1:N) {
    if (!is.nan(freq[i])){
      t <- matrix(as.integer(c(obs.1[i], obs.2[i], exp.1[i], exp.2[i])), nrow=2)
      test <- fisher.test(t)
      pvals[i] <- test$p.value
    }  
    if (progress) {
      setTxtProgressBar(pb=pb, value=i)
    }
  }
  
  result <- new("pac.result", 
                chr = data@gtdata@chromosome,
                map = data@gtdata@map,
                p.values = pvals)
  return(result)
}