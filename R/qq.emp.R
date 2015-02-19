##' @title QQ plot with empirical p-values
##' 
##' @description Plots qq plot with empirical values and empirical confidence intervals
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param obs observed values from association test 
##' @param step computation coarseness step
##' @param results vector of empirical p-values
##' @param legend if the legend is to be plotted
##' @param plot.emp plotting empirical values
##' @param conf.level confidence level (0.95 by default)
##' @return NULL
##' @examples
##'  \dontrun{
##'  qq.emp(obs=data.mm[,'Pc1df], conf=c(0.05, 0.95), result = result[,'Pc1df'],step=10, legend=T, plot.emp=T, conf.level=.95)
##'  }
##' @keywords permutations emp
##' @export qq.emp


qq.emp <- function(obs, conf=c(0.05, 0.95), result, step=10, legend=T, plot.emp=T, conf.level=.95) {
  N <- length(obs)
  obs <- -log10(obs)
  exp <- -log10(1:N/N)
  exp.s <- sort(exp)
  obs.s <- sort(obs)
  indices <- seq(from = 1, to = N, by = step)
  
  # Plotting empty graph
  plot(exp, obs, type='n', xaxt='n', yaxt='n', xlab="", ylab="", bty='n')
  xr <- par("xaxp")[1:2]
  yr <- par("yaxp")[1:2]
  abline(v=seq(xr[1], xr[2], 0.5), lty=3, col="grey")
  abline(h=seq(yr[1], yr[2], 0.5), lty=3, col="grey")
  axis(1, cex.axis=.8, col="darkgrey")
  axis(2, cex.axis=.8, col="darkgrey", las=1)
  mtext(expression(expected~~-log[10](p)), 1, 2, cex=.8, col="black") 
  mtext(expression(-log[10](p)), 2, 2, cex=.8, col="black") 
  
  # Plot empirical confidence intervals
  if (!is.null(result) & plot.emp) {
    mean_p <- apply(result, 2, mean)
    quant <- apply(result, 2, quantile, probs=conf)
    points(exp.s[indices], quant[2,indices], pch=19, cex=.5, type='l', col="black", lty=1)
    points(exp.s[indices], quant[1,indices], pch=19, cex=.5, type='l', col="black", lty=1)
    points(exp.s[indices], mean_p[indices], pch=19, cex=.5, type='l', col="darkgrey", lty=5, lwd=3)
  }
  
  # Theor p-val conf.intervals
  if (!is.null(conf)) {
    cUpp <- vector()
    cLow <- vector()
    for (i in indices) {
      cUpp[i] <- qbeta(conf[1], i, N-i+1)
      cLow[i] <- qbeta(conf[2], i, N-i+1)
    }
    points(exp[indices], -log10(cUpp[indices]), type="l", col="tomato", lty=2)
    points(exp[indices], -log10(cLow[indices]), type="l", col="tomato", lty=2)
  }
  # Theor. 
  abline(a=0, b=1, col="tomato", lty=1, lwd=1)
  # Plot observed p-values
  points(exp.s, obs.s, col="slateblue", cex=.8)
  if (!is.null(result) & !is.null(conf.level)) {
    mins <- apply(result, 1, min)
    ttest <- t.test(mins, mu=mean(mins), conf.level=conf.level)
    threshold <- ttest$conf.int[2]
    abline(h=-log10(threshold), col="tomato", lty=1, lwd=1)
    text(x = 0, -log10(threshold)+.1, labels = as.character(round(-log10(threshold), digits = 2)), cex = .8)
  }
  
  if (legend) {
    legend(x=3, y=2, 
           legend=c("obs. distr.", "emp. distr.", "theor. conf. int.", "emp. conf. int."), 
           pch=c(1,NA,NA,NA), col=c("slateblue","darkgrey", "tomato", "black"),
           lty=c(NA,5,2,1), box.lwd = 0, cex=.8)
  }
}