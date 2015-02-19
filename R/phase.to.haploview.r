##' @title Converts PHASE output to HaploView readable format.
##' 
##' @description Converts PHASE output to HaploView-friendly format.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param input a *.out file containing output from PHASE,
##' @param output filename for the result file (optional).
##' @return NULL
##' @examples
##'  \dontrun{
##'  con <- "~/Research/Behavior/results/2014-02-11_mh_aggr_chr36.phase.out"
##'  phase2haploview(input=con)
##' }
##' @keywords PHASE, HaploView, convert
##' @export phase.to.haploview

phase.to.haploview <- function(input, output=NULL) {
  require(stringr)
  f <- readLines(input)
  start <- grep("BEGIN LIST_SUMMARY",f)
  end <- grep("END LIST_SUMMARY",f)
  dat <- read.table(text = f[seq(start+1,end-1)]) 
  names(dat) <- c("haplo_id", "haplo", "count")
  tmp <- dat[,2]
  tmp <- str_replace_all(tmp, "A", "1")
  tmp <- str_replace_all(tmp, "C", "2")
  tmp <- str_replace_all(tmp, "T", "3")
  tmp <- str_replace_all(tmp, "G", "4")
  tmp <- str_replace_all(tmp, "?", "0")
  tmp <- str_replace_all(tmp, "[\\(\\)\\[\\]]", "")
  splitmp <- strsplit(tmp, "")
  splitmp <- unlist(lapply(splitmp, paste, collapse="\t"))
  dat[,2] <- splitmp

  start <- grep("BEGIN BESTPAIRS_SUMMARY",f)
  end <- grep("END BESTPAIRS_SUMMARY",f)
  dat2 <- read.table(text = f[seq(start+1,end-1)],sep=":") 
  names(dat2) <- c("id", "haplos")
  dat2$haplos <- str_replace_all(dat2$haplos, "[\\(\\)]", "")
  dat2$haplos <- str_replace_all(dat2$haplos, " ", "")
  tmp <- do.call('rbind', strsplit(dat2$haplos,",", fixed=T))
  dat2 <- cbind(dat2, hap1=as.numeric(tmp[,1]))
  dat2 <- cbind(dat2, hap2=as.numeric(tmp[,2]))
  dat2 <- dat2[,-2]
  dat2 <- cbind(dat2, hap1.seq=dat$haplo[dat2$hap1])
  dat2 <- cbind(dat2, hap2.seq=dat$haplo[dat2$hap2])
  
  if (is.null(output)) {
    output <- paste(input, ".haps", sep="")
  }
  write.table(cbind(fam.id=paste("f_", dat2$id, sep=""), 
                  id=as.character(dat2$id), 
                  hap1=as.character(dat2$hap1.seq)), 
                  file=output, row.names=F, quote=F, sep="\t")
}
phase2haploview <- phase.to.haploview
PHASE.to.Haploview <- phase.to.haploview
PHASE2Haploview <- phase.to.haploview
