##' @title Plotting a QQ-plot with desired confidence intervals. 
##' 
##' @description Plots a QQ-plot with -log10(p-value) on both axes 
##' and displays computed confidence intervals. 
##' @details Calculation of confidence intervals: 
##' The j-th order statistic from a uniform(0,1) sample has a beta(j,n-j+1) distribution. 
##' For details see Casella & Berger.
##' 
##' step goes from 1 up. Step=1 - slow and fine determination of confidence interval. 
##' Step=10 - fast and coarse determination of the confidence interval.
##' 
##' @references Casella & Berger, 2002, 2nd edition, pg 230, Duxbury.
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param pvals a vector of empirical p-values
##' @param conf a vector of the two thresholds for confidence interval 
##' @param conf.type type of confidence interval plotting, conf.type is either "shade" or "line"
##' @param conf.col color for plotting confidence interval
##' @param step computation coarseness step
##' @return NULL
##' @examples
##'  \dontrun{
##'  an0 <- qtscore(bt~antibodyLevel,data.orig)
##'  pvals <- an0[,"P1df"]
##'  plot.qq(pvals, conf=c(0.05, 0.95), conf.type="lines", conf.col="tomato", step=10)
##'  }
##' @keywords bed, gwas, p-values
##' @export plot.qq
plot.qq <- function(pvals, conf=c(0.05, 0.95), conf.type="shade", conf.col=rgb(.5,.5,1,.5), step=1) {
  lower <- conf[1]
  upper <- conf[2]
  N <- length(pvals)
  pvals <- -log10(pvals)
  expected <- -log10(1:N/N)
  maximal<- pmax(pvals, expected)
  cUpp <- vector()
  cLow <- vector()
  indices <- seq(from = 1, to = N, by = step)
  for (i in indices) {
    cUpp[i] <- qbeta(upper,i,N-i+1)
    cLow[i] <- qbeta(lower,i,N-i+1)
  }
  # QQ-plot
  qqplot(expected, pvals, ylim=c(0,max(pvals)), xlim=c(0,max(expected)), pch=19, cex=0.25, bty='n',
         las = 1, 
         xlab= expression(Expected~~-log[10](p-value)), 
         ylab=expression(Observed~~-log[10](p-value))
  ) 
  grid()
  if (conf.type == "lines") {
    points(expected[indices], -log10(cUpp[indices]), type="l", col=conf.col)
    points(expected[indices], -log10(cLow[indices]), type="l", col=conf.col)
  }
  else if (conf.type == "shade") {
    polygon(c(min(expected[indices]),max(expected[indices]),expected[indices]), 
            c(min(-log10(cUpp[indices])),max(-log10(cUpp[indices])),-log10(cLow[indices])),
            col=conf.col, border=NA)
  }
  abline(0, 1, col= 'red')
  
  res <- lm(sort(pvals)~sort(expected))
  y  <- res$coefficients[2]*expected + res$coefficients[1]  #equation of regression line 
    
 val <- pvals / y
 # val > 1 - above the line
 # val < 1 - under the line
 
 p <- c()
 for(i in 1:length(pvals)){
   if(val[i] > 10e4){
     p  <- c(p,pvals[i])
   }
 }
  
  #abline(res$coefficients[1], b = res$coefficients[2])
  mtext(text = "QQ-plot", side = 3, font = 2, cex = 1.25, line = 1.25)
}
