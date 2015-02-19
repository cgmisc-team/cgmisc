# Haplogram
setwd("~/Research/Behavior/")
data.MH1.orig <- load.gwaa.data(phenofile="data/MH_2013_canFam2.csv",
                                genofile="data/MH_2013_canFam2.raw",
                                makemap=F)
# Select only German shepherds from the MH datased
gsd <- which(data.MH1.orig@phdata$breed == "gsd" & !is.na(data.MH1.orig@phdata$curious))
data.MH1.orig <- data.MH1.orig[gsd,]
rm(gsd)
phenodata <- data.MH1.orig@phdata
rm(data.MH1.orig)
chase <- as.data.frame(phenodata$chase)
rownames(chase) <- rownames(phenodata)

setwd("~/Research/Behavior/phase/chase_phase/")
library("seqinr")
library(ape)

sequences <- read.dna(file="chase23.haps.fasta", format="fasta")
threshold <- 10
labels <- unlist(strsplit(labels(sequences), split=" "))
labels <- data.frame(name=labels[seq(1,length(labels)-1, by=2)], count=as.numeric(labels[seq(2,length(labels), by=2)]))
frequent <- which(labels$count >= threshold)
sequences <- sequences[frequent,]
d <- dist.dna(x=sequences, model="K81",as.matrix=T)
hc <- hclust(dist(d))
labels2 = as.character(labels[frequent, "name"])
hc$labels <- labels2
n = length(labels2)

op = par(bg="#ffffff", mar = c(2, 2, 2, 2))

# vector of colors
mypal = c("#556270", "#4ECDC4", "#1B676B", "#FF6B6B", "#C44D58")
N <- 4
# cutting dendrogram in N clusters
clus.cut = cutree(hc, N)
# plot
haplos <- get.bestpairs(input="chase23.out", trait=chase)

opar <- par()
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(1,1))
par(mar=c(2,0,1,0))
p <- plot(as.phylo(hc), type="phylogram", tip.color = c("olivedrab","slateblue"), 
     label.offset = 0, cex = .9, col = "red")
lastPP <- get("last_plot.phylo", envir=.PlotPhyloEnv)
#points(rep(max(lastPP$x.lim[2]),5),c(1:5), pch=19, cex=2)
plot(1,1,type='n', c(1, length(hc$order)), xlim=c(min(haplos[,"trait"],na.rm=T),max(haplos[,"trait"], na.rm=T)), yaxt='n', xaxt='n', bty='n')
axis(1, cex.axis=.7, padj=-1.5, col.axis="darkgrey", col="darkgrey")
abline(v=0, lty=2, col="grey")
tmp <- clus.cut[hc$order]
for (i in unique(tmp)) {
  #data <- rnorm(1000, sample(seq(-.3,.3,length.out=100),size=1), 1)
  haplo <- names(tmp[tmp==i])
  data <- haplos[which(haplos$hap1 %in% haplo | haplos$hap1 %in% haplo),"trait"]
  y <- which(tmp %in% i)[1]
  if (y %% 2 != 0) {
    colour <- "olivedrab"
  } else {
    colour <- "slateblue"
  }
  my.boxplot(data, y=y, cols=c(colour, "grey"))
}

my.boxplot <- function(data, y, cols=c("slateblue","grey")) {
  q <- quantile(data, probs=c(.25, .5, .75))
  t1 <- q[1] - 1.5 * IQR(data)
  t2 <- q[3] + 1.5 * IQR(data)
  outl.t1 <- data[data < t1]
  outl.t2 <- data[data > t2]
  lines(x=c(min(data[data >= t1]), max(data[data <= t2])), y=c(y,y), col="darkgrey")
  points(mean(data), y, pch=3, col="black", cex=2)
  lines(x=c(q[1], q[3]), y=c(y, y), col=cols[1], lwd=10)
  points(median(data), y, pch=21, bg=cols[2])
  points(outl.t1, rep(y, length(outl.t1)), pch=1, cex=.7, col=cols[1])
  points(outl.t2, rep(y, length(outl.t2)), pch=1, cex=.7, col=cols[1])
}

get.bestpairs <- function(input, trait) {
  require(stringr)
  f <- readLines(input)  
  start <- grep("BEGIN BESTPAIRS_SUMMARY", f)
  end <- grep("END BESTPAIRS_SUMMARY", f)
  dat2 <- read.table(text = f[seq(start+1,end-1)],sep=":") 
  names(dat2) <- c("id", "haplos")
  dat2$haplos <- str_replace_all(dat2$haplos, "[\\(\\)]", "")
  dat2$haplos <- str_replace_all(dat2$haplos, " ", "")
  tmp <- do.call('rbind', strsplit(dat2$haplos,",", fixed=T))
  dat2 <- cbind(dat2, hap1=paste("haplo", tmp[,1], sep=""))
  dat2 <- cbind(dat2, hap2=paste("haplo", tmp[,2], sep=""))
  dat2 <- dat2[,-2]
  tmp <- trait[dat2$id,]
  dat2 <- cbind(dat2, trait=tmp)
  dat2
}

