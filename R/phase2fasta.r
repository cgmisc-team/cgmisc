##' @title Converts haplotypes from PHASE output to FASTA format.
##' 
##' @description Returns a FASTA file with all PHASE-derived haplotypes.
##' @author Marcin Kierczak \email{Marcin.Kierczak@@imbim.uu.se}
##' @param input output from PHASE (.out file),
##' @param output output filename (optional).
##' @return NULL
##' @examples 
##' \dontrun{
##'  con <- "~/Research/Behavior/results/2014-02-11_mh_aggr_chr36.phase.out"
##'  phase2fasta(input=con)
##' }
##' @keywords fasta, gwas, p-values
##' @export phase2fasta
phase2fasta <- function(input, output=NULL, filter=NULL) {
  require(stringr)
  f <- readLines(input)
  start <- grep("BEGIN LIST_SUMMARY",f)
  end <- grep("END LIST_SUMMARY",f)
  dat <- read.table(text = f[seq(start+1,end-1)]) 
  names(dat) <- c("haplo_id", "haplo", "count")
  if (!is.null(filter)) {
    dat <- dat[dat$count >= filter, ]
  }
  headline <- paste(">haplo_", dat$haplo_id, "_count_", dat$count, sep="")
  result <- paste(headline,"\n",dat$haplo, sep="")
  if (is.null(output)) {
    output <- paste(input, ".fasta", sep="")
  }
  writeLines(result, con=output, sep="\n")
}