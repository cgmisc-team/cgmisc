##' @title QQ plot with empirical p-values
##' 
##' @description Plots qq plot with empirical values and empirical confidence intervals
##' @author Marcin Kierczak <\email{Marcin.Kierczak@@imbim.uu.se}>
##' @param data a gwaa.data-class object used to fit the original model
##' @param obs observed values from association test 
##' @param emp vector of empirical p-values or  a list, result of running polygenic model
##' @param N if polygenic model supplied, a number of permutations to run
##' @param step computation coarseness step
##' @param legend if the legend is to be plotted
##' @param plot.emp plotting empirical values
##' @param conf.level confidence level (0.95 by default)
##' @param show.pb a logical indicating whether progress bar will be shown
##' @conf a vector defininy the lower and the upper confidence interval, default 5% CI.
##' @return NULL
##' @examples
##'  \dontrun{
##'  qq.emp(obs=data.mm[,'Pc1df], conf=c(0.025, 0.975), emp = result[,'Pc1df'], step=10, legend=T, plot.emp=T, conf.level=.95)
##'  }
##' @keywords permutations emp
##' @export qq.emp
plot.qq <- function(data=NULL, obs, emp, N=30, step=10, legend=T, plot.emp=T, conf.level=.95, conf=c(0.025, 0.975), show.pb=T) {
  require(GenABEL)
  N <- length(obs)
  obs <- -log10(obs)
  exp <- -log10(1:N/N)
  exp.s <- sort(exp)
  obs.s <- sort(obs)
  indices <- seq(from = 1, to = N, by = step)
  
  # Plotting empty graph
  plot(exp, obs, type='n', xaxt='n', yaxt='n', xlab="", ylab=expression(observed~~-log[10](p)), 
       cex.lab=.8, bty='n', mgp=c(2, 1, 0))
  xr <- par("xaxp")[1:2]
  yr <- par("yaxp")[1:2]
  # abline(v=seq(xr[1], xr[2], 0.5), lty=3, col="grey")
  # abline(h=seq(yr[1], yr[2], 0.5), lty=3, col="grey")
  grid()
  axis(1, cex.axis=.8, col="darkgrey")
  axis(2, cex.axis=.8, col="darkgrey", las=1)
  mtext(expression(expected~~-log[10](p)), 1, 2, cex=.8, col="black") 
  # mtext(expression(-log[10](p)), 2, 2, cex=.8, col="black") 
   
  # Perform permutations if necessary
  if (class(emp) == "polygenic") {
    if (is.null(data)) {
      stop("Error! gwaa.data not supplied!")
    }
    perm.result <- c()
    h2h <- emp
    h2h.tmp <- emp
    if (show.pb) {
      pb <- txtProgressBar(min=0, max=N, initial=0, style=3)
    }
    for (i in 1:N) {
      h2h.tmp$grresidualY <- sample(h2h$grresidualY)
      tmp <- qtscore(h2h.tmp$grresidualY, data, clambda = F)
      #perm.result <- rbind(perm.result, sort(-log10(tmp@results$P1df)))
      if (show.pb) { 
        setTxtProgressBar(pb, i)
      }
    }
  emp <- perm.result
  }
  
  # Plot empirical confidence intervals
  if (!is.null(emp) & plot.emp) {
    mean_p <- apply(emp, 2, mean)
    quant <- apply(emp, 2, quantile, probs=conf)
    points(exp.s[indices], quant[2,indices], pch=19, cex=.5, type='l', col="black", lty=1)
    points(exp.s[indices], quant[1,indices], pch=19, cex=.5, type='l', col="black", lty=1)
    points(exp.s[indices], mean_p[indices], pch=19, cex=.5, type='l', col="black", lty=5, lwd=1)
  }
  
  # Theor p-val conf.intervals
  if (!is.null(conf)) {
    cUpp <- vector()
    cLow <- vector()
    for (i in indices) {
      cUpp[i] <- qbeta(conf[1], i, N-i+1)
      cLow[i] <- qbeta(conf[2], i, N-i+1)
    }
    points(exp[indices], -log10(cUpp[indices]), type="l", col="tomato", lty=1)
    points(exp[indices], -log10(cLow[indices]), type="l", col="tomato", lty=1)
  }
  
  # Theoretical distribution 
  abline(a=0, b=1, col="tomato", lty=5, lwd=1)
  
  # Plot observed p-values
  # points(exp.s, obs.s, col="slateblue", cex=.7, type='l')
  points(exp.s, obs.s, col="slateblue", cex=.5, pch=19)
  
  if (!is.null(emp) & !is.null(conf.level)) {
    mins <- apply(emp, 1, min)
    ttest <- t.test(mins, mu=mean(mins), conf.level=conf.level)
    threshold <- ttest$conf.int[2]
    abline(h=-log10(threshold), col="darkgrey", lty=1, lwd=1)
    mtext(text = round(-log10(threshold), digits = 1), 2, line = 0, cex = .7, adj = 1.5, , at = -log10(threshold), col="darkgrey")
  }
  
  if (legend) {
    legend("topleft", 
           legend=c("theor. distr.","theor. conf. int.", "emp. distr.","emp. conf. int.", "emp. thr."), 
           pch=c(NA,NA,NA,NA,NA), col=c("tomato","tomato","black","black", "darkgrey"),
           lty=c(5,1,5,1,1), box.lwd = 0, cex=.8)
  }
}

