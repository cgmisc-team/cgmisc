setwd("~/Research/Behavior/")
load("data/varexp.data.Rd")
data <- data.tmp
rm(data.tmp)
library(cgmisc)
library(wesanderson)

chr.list <- c(1,2,3)
pop <- sample(c(1, 2), size=nids(data), replace=T)
x <- c()
for (i in 1:length(chr.list)) {
  chr <- data[,data@gtdata@chromosome == chr.list[i]]
  pacs <- pop.allele.counts(data=chr, pops=pop, progress=T)
  #W <- get.windows(chr, d=10^6)
  W <- get.overlapping.windows(data = chr, w=125e3, o=25e3)[[2]]
  y <- -log10(pacs@p.values)
  means <- get.window.means(y, W)
  x <- c(x, means)
}

dat <- data[,data@gtdata@chromosome %in% chr.list]
result <- data.frame(mean.p.value=x, chr=dat@gtdata@chromosome, map=dat@gtdata@map)
result.sorted <- result[with(result, rev(order(mean.p.value))), ]

# Plot averaged results
col <- wes.palette(n=2, name="Royal1")[1 + as.numeric(dat@gtdata@chromosome) %% 2]
plot(x, type='h', bty='n', col=col, las=1, xaxt='n', ylab="mean -log10(p-value)", xlab="chromosome")
midpoints <- get.chr.midpoints(data=dat)
axis(side=1, at=midpoints, labels=chr.list)

# Determine the confidence interval and plot it as a line
#conf.lvl <- 0.95 # Set desired conf. int. if you give 1%, threshold will be plotted at 99%.
#t <- t.test(x, mu=mean(x), conf.level=conf.lvl)
#abline(h=t$conf.int[2], lty=3)

# Threshold based on 1% top values
q99 <- quantile(x=x, probs=c(0.99))
abline(h=q99, lty=3)

# Zooming-in
head(result.sorted, n=10) # Top 10 results
region <- c(1, 911277)
chr <- 1
zoom <- which(data@gtdata@chromosome == chr & data@gtdata@map >= region[1] & data@gtdata@map <= region[2])
# You can change type to 'h' for histogram-like plot
plot(result[zoom, "map"], x[zoom], ylim=c(0, max(x)), xlab="Mbp", ylab="mean -log10(p-value)", bty='n', las=1, xaxt='n', type='n')
grid()
points(result[zoom, "map"], x[zoom], col="slateblue", type='l')
#abline(h=t$conf.int[2], lty=3)
abline(h=q99, lty=2)
axis(side=1, at=pretty(result[zoom, "map"]), labels=pretty(result[zoom, "map"])/10e6)

get.chr.midpoints <- function(data) {
  chromosomes <- unique(as.numeric(data@gtdata@chromosome))
  if (length(chromosomes) > 1) {
    chr.mid <- c()
    for (chr in chromosomes) {
      coord <- which(as.numeric(data@gtdata@chromosome) == chr)
      chr.mid <- c(chr.mid, floor((coord[1] + coord[length(coord)])/2))
    }
  }
return(chr.mid)
}

makemap <- function(data) {
  a <- data@gtdata
  cat("increase in map order FORCED\n")
  chun <- levels(a@chromosome)
  if (any(chun != "X")) {
    numchun <- sort(as.numeric(chun[chun != "X"]))
    gsize <- max(a@map[a@chromosome == as.character(numchun[1])]) / 5
    if (length(numchun) > 1) {
      for (i in c(2:(length(numchun)))) {
        inc <- max(a@map[a@chromosome == as.character(numchun[i - 1])]) + gsize
        a@map[a@chromosome == as.character(numchun[i])] <- a@map[a@chromosome == as.character(numchun[i])] + inc
      }
    }
    if (any(chun == "X")) {
      inc <- max(a@map[a@chromosome == as.character(numchun[length(numchun)])]) + gsize
      a@map[a@chromosome == "X"] <- a@map[a@chromosome == "X"] + inc
    }
  }
  return(a@map)
}

get.window.means <- function(x, W) {
  X <- matrix(rep(x, times=length(x)), nrow=length(x), byrow=T)
  X[W == F] <- NA
  means <- apply(X, 1, mean, na.rm=T)
  return(means)
}

get.windows <- function(data, d=10^6) {  
  m <- as.matrix(data@gtdata@map)
  dm <- as.matrix(dist(m, diag=T))
  dw <- matrix(F, nrow=dim(dm)[1], ncol=dim(dm)[1])
  dw[dm <= d] <- T
  return(dw)
}

dat <- data[,data@gtdata@chromosome == 1]

get.overlapping.windows <- function(data, w = 125e3, o = 25e3) {
  hw <- floor(0.5 * w)  
  m <- as.matrix(data@gtdata@map)
  map <- data@gtdata@map
  centers <- seq(map[1] + hw, rev(map)[1] - hw, by=(w-o))
  starts <- centers - hw
  stops <- centers + hw
  windows <- data.frame(start=starts, midpoint=centers, stop=stops)
  M <- matrix(FALSE, nrow = length(starts), ncol=length(dat@gtdata@map))
  for (i in 1:length(starts)) {
    M[i, (dat@gtdata@map >= starts[i] & dat@gtdata@map <= stops[i])] <- T
  }
  return(list(windows, M))
}

tmp <- get.overlapping.windows(dat)

