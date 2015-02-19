setwd("~/Dropbox/20140217/")
library("GenABEL")
library("compiler")
require("fields")
enableJIT(level=3)
data1 <- load.gwaa.data(phenofile="AtopyGWAS_20120628_pool3_b.pheno", genofile="20120628-atopy.raw")

tmp <- data1[,data1@gtdata@chromosome=="27"]
allowed.na <- .2
allowed.het <- .1
pb <- txtProgressBar(min=0, max=nsnps(tmp), initial=0, style=3)
result <- data.frame(matrix(0, nrow=nids(tmp), ncol=nsnps(tmp)), row.names=tmp@gtdata@idnames)
colnames(result) <- tmp@gtdata@snpnames
for (i in 1:nsnps(tmp)) {
  window <- get.adjacent.markers(tmp, marker=i, size.bp=25e3)
  roh.test <- apply(X=window, 1, is.roh, allowed.na, allowed.het)
  inds <- names(roh.test[roh.test==1])
  if (length(inds) > 0) {
    snps <- colnames(window)
    result[inds, snps] <- result[inds, snps] + 1
  }
  setTxtProgressBar(pb,i)
}
res <- result/rowSums(result)
res <- res > 0.05

is.roh <- function(x, allowed.na, allowed.het) {
    if (allowed.na < 1) {
      allowed.na <- floor(allowed.na * length(x))
    }
    if (allowed.het < 1) {
      allowed.het <- floor(allowed.het * length(x))
    } 
    is.roh <- 0
    tmp <- x == 1
    if (sum(is.na(x)) <= allowed.na & sum(x, na.rm=T) <= allowed.het) {
      is.roh <- 1
    }
    return(is.roh)
}













hr <- het.scan.ind(data1, 1, 27)
het.scan.ind <- function(data, ind=1, chr, size=25e3, het.threshold=1) {
  subset <- data[ind, data@gtdata@chromosome==chr] 
  scan <- data.frame(map=subset@gtdata@map, 
                     windows=rep(0, 
                     times=nsnps(subset)), 
                     het=rep(0, times=nsnps(subset)))
  rownames(scan) <- subset@gtdata@snpnames
  pb <- txtProgressBar(min=0, max=nsnps(subset), initial=0, style=3)
  for (i in 1:nsnps(subset)) {
    tmp <- get.adjacent.markers(subset, subset@gtdata@snpnames[i], size)
    hetvec <- sum(summary(tmp)$Q.2 > 0)
    if (hetvec <= het.threshold) {
      scan[tmp@snpnames,3] <- scan[tmp@snpnames,3] + 1
    }
    scan[tmp@snpnames,2] <- scan[tmp@snpnames,2] + 1
    setTxtProgressBar(pb,i)
  }
  prop <- scan$het/sum(scan$het)
  scan <- cbind(scan, prop)
  return(scan)
}
plot(scan$map, scan$het/scan$windows, pch=19, cex=.5)

get.adjacent.markers <- function(data, marker, size.bp=1e3) {
  chr <- data@gtdata@chromosome[marker]
  chr.data <- data[,which(data@gtdata@chromosome == chr)]
  coord <- chr.data@gtdata@map[marker]
  start <- pmax(min(chr.data@gtdata@map), coord - size.bp)
  end <- pmin(coord + size.bp, max(chr.data@gtdata@map))
  tmp <- chr.data[,chr.data@gtdata@map >= start & chr.data@gtdata@map <= end]
  return(tmp@gtdata)
}
