##' @title Compute Fst fixation index for two populations
##' 
##' @description Fixation index Fst is a measure of population differentiation due to genetic structure. 
##' Given a set of genotypes in two populations, the function computes fixation index 
##' Fst and corresponding indices: Fit and Fis.
##' @param data a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}}
##' @param pops a vector of two values indicating to which population an individual belongs
##' 
##' @details 
##' \itemize{
##'   \item data -- a standard \code{\link[GenABEL]{gwaa.data-class}} object
##'   \item pops -- a vector of two values indicating to which population an individual belongs. 
##'   Typically, one uses a vector of zeroes and ones where 0 marks an individual belonging to 
##'   population 1 and 1 marks an individual belonging to population 2. Often, the vector is a result 
##'   of clustering in MDS-scaled genomic kinship space
##' }
##' @return an \code{\link[cgmisc]{fstats.result}} class object
##' @references Holsinger, Kent E.; Bruce S. Weir (2009). "Genetics in geographically structured populations: 
##' defining, estimating and interpreting FST". Nat Rev Genet 10 (9): 639-650. doi:10.1038/nrg2611. ISSN 1471-0056. PMID 19687804.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@slu.se}>
##' @examples 
##' \dontrun{
##' fstats  <- compute.Fstats(data, pops)
##' }
##' @seealso \code{\link[GenABEL]{gwaa.data-class}}, \code{\link[cgmisc]{fstats.result}}
##' @keywords FST fixation index population structure heterozygosity
##' @export compute.fstats
##' 
compute.fstats <- function(data, pops) {
  # Copyright Marcin.Kierczak@slu.se
  # v. 2.0, 17.01.2013
  pop.names <- unique(pops)
  require("GenABEL")
  geno.double <- as.double.snp.data(data@gtdata)
  subpop_data <- list()
  sum.obs.het <- rep(0, times=dim(geno.double)[2])
  sum.exp.het <- rep(0, times=dim(geno.double)[2])
  sum.N <- rep(0, times=dim(geno.double)[2])
 
  for (pop_nr in c("global", pop.names)) {
    if (pop_nr == "global") {
      tmp_data <- geno.double
    }
    else {
      tmp_data <- geno.double[pops == pop_nr,]
    }
    sum1 <- colSums(tmp_data == 0, na.rm = T)
    sum2 <- colSums(tmp_data == 1, na.rm = T)
    sum3 <- colSums(tmp_data == 2, na.rm = T)
    subpop_data[[pop_nr]] <- data.frame(AA = sum1, Aa = sum2, aa = sum3)
    tmp <- subpop_data[[pop_nr]]
    tmp[,"N"] <- colSums(t(tmp))
    tmp[,"p"] <- (2 * tmp["AA"] + tmp["Aa"]) / (2 * tmp["N"])
    tmp[,"q"] <- 1 - tmp[,"p"]
    tmp[,"exp.AA"] <- tmp[,"N"] * tmp[,"p"]^2 
    tmp[,"exp.Aa"] <- tmp[,"N"] * 2 * tmp[,"p"] * tmp[,"q"]
    tmp[,"exp.aa"] <- tmp[,"N"] * tmp[,"q"]^2
    tmp[,"obs.het"] <- tmp[,"Aa"] / tmp[,"N"]
    tmp[,"exp.het"] <- 2 * tmp[,"p"] * tmp[,"q"]
    tmp[,"dev.het"] <- tmp[,"exp.het"] - tmp[,"obs.het"]
    # Inbreeding coeff. Fs
    tmp[,"Fs"] <- tmp[,"dev.het"] / tmp[,"exp.het"]
    subpop_data[[pop_nr]] <- tmp
    sum.obs.het <- sum.obs.het + tmp[,"obs.het"] * tmp[,"N"]
    sum.exp.het <- sum.exp.het + tmp[,"exp.het"] * tmp[,"N"]
    sum.N <- sum.N + tmp[,"N"]
  }
  # p-bar and q-bar (freq of allele A/a over the total pop)
  p.bar <- subpop_data[["global"]][,"p"] 
  q.bar <- subpop_data[["global"]][,"q"]
  # Global heterozygosity indices
  HI <- sum.obs.het/sum.N
  HS <- sum.exp.het/sum.N
  HT <- 2 * p.bar * q.bar
  # Global F-statistics
  FIS <- (HS - HI) / HS 
  FST <- (HT - HS) / HT
  FIT <- (HT - HI) / HT
  glob <- data.frame(p.bar=p.bar, q.bar=q.bar, HI=HI, HS=HS, HT=HT, FIS=FIS, FST=FST, FIT=FIT)
  result <- new("fstats.result", pops=subpop_data, glob=glob)
  result
}
compute.Fstats <- compute.fstats
