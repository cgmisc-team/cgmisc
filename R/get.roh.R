##'@title Runs of homozygosity
##'@description Identifies runs of homozygosity in overlapping windows across chromosome.
##'@author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>, Jagoda Jablonska <\email{jagoda100jablonska@@gmail.com}>
##'@param data a gwaa.data class object as used by \code{\link[GenABEL]{gwaa.data-class}}
##'@param chr number of chromosome 
##'@param LW a matrix of windows coordinares returned by \code{get.overlapping.windows} function
##'@param hetero.zyg a matrix of average heterozygosities for windows returned by \code{het.for.overlap.wind} function
##'@param threshold threshold of homozygosity. All windows with homozygosity below this limit will be treated as homozygous.
##'@param strict if one window with lower homozygosity should be tolerated
##'@return matrix with columns respectively : 1) number of first window of run, 2) coordinate of first window in a run,
##'3) coordinate of last window in a run, 4) length of run (number of windows)
##'@seealso \code{\link[cgmisc]{get.overlapping.windows}} \code{\link[cgmisc]{het.overlap.wind}}
##'@keywords heterozygosity, homozygosity, overlapping window 
##'
##'@export get.roh

get.roh  <- function (data, chr, LW, hetero.zyg, threshold = 0.23, strict = TRUE) {
  ## calculating threshold for ROH
  gtype = data@gtdata
  het.all  <- descriptives.marker(data)[[7]]
  n = gtype@nsnps - gtype@nids
  L = log(0.05/n) / log(1-het.all) #minimal number of SNPs
  hetero.zyg <- na.omit(hetero.zyg)
  
  ## identifying regions with high homozygosity
  temp.roh  <- list()
  i = 1
  while (i < nrow(hetero.zyg)) {
    if(i >= nrow(hetero.zyg)){break}
    if(hetero.zyg[i,2] <= threshold) {
      het.prev  <- hetero.zyg[i,2]
      temp.roh[i]  <- TRUE
      for (y in i+1:nrow(hetero.zyg)){
        if (hetero.zyg[y,2] <= threshold){
          temp.roh[y]  <- TRUE 
          het.prev  <- hetero.zyg[y,2]
          if (y >= nrow(hetero.zyg)) {
            temp.roh -> final
            i = i+1
            break 
          } else {
            next
          }
        } else {
          temp.roh[y]  <- FALSE
          i  = y+1
          break 
        }
      }
    } else {
      temp.roh[i]  <- FALSE
      i = i+1
      next 
    }
  }
  temp.roh  -> final
  
  ## Identifying runs of homozygosity across the windows
  dist  <- as.matrix(as.data.frame(LW[1]))
  logical  <- unlist(final)
  if(strict == FALSE){
  for (k in 2:length(logical))
  {
    if (logical[k] == FALSE & logical[k+1] == TRUE & logical[k-1] == TRUE) {
      logical[k] <- TRUE
    }
  }
  }
  blocks  <- c(logical[1] == TRUE, diff(logical))
  runs  <- matrix(ncol = 4, nrow = 5000)
  colnms  <- c('window','begin','end','length')
  colnames(runs) <- paste(colnms)
  index = 0
  i = 1
  while( i < length(blocks)) {
    if( blocks[i] == 1){
      index = index +1
      runs[index,1]  <- i
      runs[index,2]  <- dist[i,1]
      while (blocks[i] != -1 & !is.na(blocks[i])){
        i = i+1
      }
      runs[index, 3]  <- dist[i-1,3]
      beg  <- runs[index,2]
      end  <- runs[index,3]
      nsnp  <- length(which(data@gtdata@chromosome == chr & data@gtdata@map > beg & data@gtdata@map < end ))
      nsnp -> runs[index,4]
      next
    } else {
      i = i+1
      next
    }     
  }
  runs  <-  runs[rowSums(is.na(runs)) != ncol(runs),]
  runs <- runs[which(runs[,'length'] > L ),]
  return (runs)
}
