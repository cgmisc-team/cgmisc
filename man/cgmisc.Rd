\name{cgmisc-package}
\alias{cgmisc-package}
\alias{cgmisc}
\docType{package}
\title{
cgmisc: Enhanced Genome-wide Association Analyses and Visualisation.
}
\description{
Package faciliating and enhancing data analyses and visualisation in genome-wide association studies. Contains functions for analyses of population structure, linkage structure, heterozygosity as well as wide range of data subsetting and manipulation.
}
\details{
\tabular{ll}{
Package: \tab cgmisc\cr
Version: \tab 2.9.9\cr
Date: \tab 2015-04-28\cr
Depends: \tab 
R (>= 3.0.1),
GenABEL (>= 1.7.4),
GenomicRanges (>= 1.18.4),
rtracklayer (>= 1.26.2),
ggplot2 (>= 1.0.1),
grid (>= 3.1.1),
wq (>= 0.4.3)\cr
VignetteBuilder: \tab knitr\cr
License: \tab GPL (>=2.10)\cr
}

Index:
\tabular{ll}{
\code{\link{boxplot.snp}} \tab Box plot with genotypes \cr
\code{\link{boxplot.snp.twoWay}} \tab SNP_box_twoWay \cr
\code{\link{choose.top.snps}} \tab Select the highest r2 SNPs to a given top marker. \cr
\code{\link{chr.x.fix.canfam}} \tab Fix canine X chromosome naming. \cr
\code{\link{clump.markers}} \tab Clump markers according to their LD. \cr
\code{\link{compute.fstats}} \tab Compute Fst fixation index for two populations. \cr
\code{\link{create.haploview.info}} \tab Prepare a PHASE input file from gwaa data. \cr
\code{\link{epistasis.scan}} \tab Scan for epistasis between a specified SNP and all others.\cr
\code{\link{get.adjacent.markers}} \tab Get adjacent markers within a distance around the given marker. \cr
\code{\link{get.chr.midpoints}} \tab Get middle point for each chromosome. \cr
\code{\link{get.erv}} \tab Get endogenous retroviral sequences in canine genome, \cr
\code{\link{get.overlap.windows}} \tab Get overlapping windows. \cr
\code{\link{get.roh}} \tab Identify runs of homozygosity. \cr
\code{\link{gwaa.to.bed}} \tab Create a BED file with GWAS p-values in a region. \cr
\code{\link{gwaa.to.bigrr}} \tab Converts gwaa-data-class GenABEL object to a bigRR-readable format. \cr
\code{\link{gwaa.to.phase}} \tab Prepare a PHASE input file from gwaa data.\cr
\code{\link{gwaa.to.vgwas}} \tab Convert GenABEL data structure to vGWAS data structure. \cr
\code{\link{het.overlap.wind}} \tab Calculates heterozygosity for windows. \cr
\code{\link{open.region.ucsc}} \tab Open UCSC GenomeBrowser for selected genomic region. \cr
\code{\link{phase.to.fasta}} \tab Converts haplotypes from PHASE output to FASTA format. \cr
\code{\link{phase.to.haploview}} \tab Converts PHASE output to HaploView readable format. \cr
\code{\link{plot.clumps}} \tab Plot marker clumps on Manhattan plot. \cr
\code{\link{plot.erv}} \tab Plot endogenous retroviral sequences in canine genome. \cr
\code{\link{plot.fstats}} \tab Plot results of 'compute.fstats'. \cr
\code{\link{plot.genes}} \tab Plot genes provided in a bed file. \cr
\code{\link{plot.ld.decay}} \tab Plots LD decay pattern for a given data. \cr
\code{\link{plot.manhattan.genes}} \tab Plot LD pattern and MAF in a Manhattan plot. \cr
\code{\link{plot.manhattan.ld}} \tab Plot LD pattern and MAF in a Manhattan plot. \cr
\code{\link{plot.overlap}} \tab Plots average heterozygosities for overlapping windows. \cr
\code{\link{plot.pac}} \tab Plot p-values of allele count differences between populations. \cr
\code{\link{plot.qq}} \tab QQ plot with empirical p-values. \cr
\code{\link{pop.allele.counts}} \tab Compare allele counts in two populations. \cr
}
}
\author{
Marcin Kierczak <marcin.kierczak@imbim.uu.se>\cr
Jagoda Jablonska <jagoda100jablonska@gmail.com>\cr
Simon Forsberg <simon.forsberg@slu.se>\cr
Veronika Scholz <veronika.scholz@imbim.uu.se>

Maintainer: Marcin Kierczak <marcin.kierczak@imbim.uu.se>
}
\references{
Soon}
\keyword{ gwas }
\seealso{

\code{\link[<pkg>:<pkg>-package]{GenABEL}}
}
