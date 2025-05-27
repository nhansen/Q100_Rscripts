setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/PaperFigures/Figures/VCFgenome")

library(colorspace)
library(Hmisc)

# no scientific notation on axes
options(scipen=5)

errorstats <- read.table("totalhcbasesandvariantcounts.txt", sep="\t")
names(errorstats) <- c("ref", "mingq", "covered", "totalerrors", "snperrors", "indelerrors")

errorstats$percentcovered <- errorstats$covered/5999408148
errorstats$totalqv <- -10.0*log10(errorstats$totalerrors/errorstats$covered)
errorstats$snpqv <- -10.0*log10(errorstats$snperrors/errorstats$covered)
errorstats$indelqv <- -10.0*log10(errorstats$indelerrors/errorstats$covered)

grch38errorstats <- errorstats[errorstats$ref=="GRCh38",]
chm13errorstats <- errorstats[errorstats$ref=="CHM13",]

qualcolors <- c("#44AA99", "#332288", "#882255", "#888888")
barcolors <- sapply(qualcolors, function(x) {c(x, x)})

barplot(c(grch38errorstats$percentcovered, chm13errorstats$percentcovered), main="Percent of benchmark covered by VCF-constructed genomes", names.arg=c("GRCh38", "CHM13"), beside=TRUE, col=t(barcolors))
barplot(ploterrorcountsbyqual(nocensattotalerrors), names.arg=c("CHM13", "GRCh38"), beside=TRUE, col=t(barcolors), main="Total errors outside censat regions for Revio calls")
legend("topright", c("10", "20", "30", "40"), col=qualcolors, pch=15, title="Min GQ", cex=1.4)

variantcounts <- read.table("totalvariantcounts.txt", sep=" ")
hcbasecounts <- read.table("totalhighconfbases.txt", sep=" ")

barplot(ploterrorcountsbyqual(variantcounts), names.arg=c("CHM13", "GRCh38"), beside=TRUE, col=t(barcolors), ylab="Variant Count", main="Variant calls with GQ over threshhold")
barplot(ploterrorcountsbyqual(hcbasecounts), names.arg=c("CHM13", "GRCh38"), beside=TRUE, col=t(barcolors), ylab="High Conf Bases", main="Bases called with GQ over threshhold")
legend("topright", c("10", "20", "30", "40"), col=qualcolors, pch=15, title="Min GQ")

# statistics re: GQ values in gVCF file (from Revio PacBio/deepvariant pipeline run)
refcallgqs <- read.table("gVCFs/GRCh38.gVCF.RefCall.GQvals.txt")
refregiongqs <- read.table("gVCFs/GRCh38.gVCF.RefRegions.GQvals.txt")
varcallgqs <- read.table("gVCFs/GRCh38.gVCF.VarCalls.GQvals.txt")

# reference regions with GQ values:
refregionswithgqs <- read.table("gVCFs/GRCh38.gVCF.RefRegions.bed")
refregiongqbases <- sapply(seq(1,50), function(gq) {sum(refregionswithgqs[refregionswithgqs$V4==gq, "V3"]-refregionswithgqs[refregionswithgqs$V4==gq, "V2"] + 1)})
