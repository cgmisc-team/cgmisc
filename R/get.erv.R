##' @title Get endogenous retroviral sequences in canine genome
##' 
##' @description Returns a list of endogenous retroviral sequences (ERV) in a predefined chunk of the canine genome.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>, Patric Jern <\email{Patric.Jern@@imbim.uu.se}>
##' @param chr chromosome, e.g. "chr17", 
##' @param coords a vector of coordinates,
##' @param src source database, currently only canFam3 is supported.
##' @return a dataframe containing the following columns: id, chromosome, strand, start, end, length, score, subgenes (names) 
##' @examples
##'  \dontrun{get.erv("chr16", coords=c(1e6,2e6))
##'  }
##' @keywords get, ERV
##' @export

get.erv <- function(chr=NA, coords=c(NA,NA), src="canFam3cgmisc.db") {
  require("RSQLite")
  src <- paste("db/", src, sep="")
  if (is.na(chr)) {
    stop("Chromosome not valid.")
  }
  sqlite    <- dbDriver("SQLite")
  db <- dbConnect(sqlite, src)
  query <- ""
  query <- paste("SELECT * FROM ERV WHERE chromosome == '", chr, "'", sep = "")
  if (!is.na(coords[1])) {
    query <- paste(query, " AND start >= ", coords[1], sep = "")
  }
  if (!is.na(coords[2])) {
    query <- paste(query, " AND end <= ", coords[2], sep = "")
  }
  result <- dbGetQuery(db, query)
  tmp <- dbDisconnect(db)
  return(result)
}
