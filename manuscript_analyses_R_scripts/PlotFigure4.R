setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/PaperFigures/Figures")

library(colorspace)
library(Hmisc)

source("/Users/nhansen/OneDrive/HG002_diploid_benchmark/Q100_Rscripts/manuscript_analyses_R_scripts/ReadBenchComparisonPlotFunctions.R")

################
### FIGURE 4 ###
### Sequencing coverage and accuracy of various platforms compared against ###'
### the HG002 genome benchmark. [Essential read quality statistics from ###
### Nanopore, PacBio, Illumina, Element, etc. illustrating subtle base ###
### call and coverage biases not evident from typically reported QV scores. ###
### Things like real vs. reported QV values, indel biases, coverage ###
### uniformity, transition/transversion biases ###
################

outdir <- "Figure4Output"

# currently plotted read sets (alternatives below)
readsetnames <- c("ont_epi2me_q28", "hifi_revio_pbmay24", "element_ultraq_jun2024", "illumina_googlepcrfree")
platformlabels <- c("ONT Q28", "HiFi Revio", "Element", "Illumina")

#readsetnames <- c("ont_epi2me_q28", "lc24_apk", "hifi_revio_pbmay24", "NIST_onso_2024Q1", "element_ultraq_jun2024", "illumina_2x250")
#readsetnames <- c("ont_epi2me_q28", "lc24_apk", "hifi_revio_pbmay24", "element_ultraq_jun2024", "illumina_2x250")
#readsetnames <- c("ont_epi2me_q28", "lc24_apk", "hifi_revio_pbmay24", "element_ultraq_jun2024", "illumina_googlepcrfree_ds0.1")
#platformlabels <- c("ONT Q28", "ONT APK", "HiFi Revio", "Onso", "Element", "Illumina")
#platformlabels <- c("ONT Q28", "ONT APK", "HiFi Revio", "Element", "Illumina")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
readplatformcolors <- c("#6699CC", "#01541F", "#332288", "#661100", "#5D5D5D")
platformlinetype <- c(1, 2, 3, 4, 5, 6)
platformpchvals <- c(0, 1, 2, 5, 6, 8)
filledplatformpchvals <- c(15, 16, 17, 23, 25, 8)
title <- ""

pdf("Figure4Output/QVAccuracy.pdf")
read_qv_plot(readsetnames, platformlabels)
dev.off()

### Make mononucleotide accuracy plot ###

pdf("Figure4Output/ReadMononucAccuracy.pdf")  
read_mononucacc_plot(readsetnames, platformlabels)
dev.off()

pdf("Figure4Output/ReadMononucErrorRate.pdf")  
read_mononucerror_plot(readsetnames, platformlabels)
dev.off()

pdf("Figure4Output/ReadMononucQVScores.pdf")  
read_mononucqvscore_plot(readsetnames, platformlabels)
dev.off()

# Plot substitution rate-by-type histogram

pdf("Figure4Output/ReadSubstitutionRates.pdf", width=11, height=6)
read_substitutions_plot(readsetnames, platformlabels, outputdir="Figure4Output", legend=TRUE)
dev.off()

### Indel error plot ###

pdf("Figure4Output/ReadIndelRates.pdf")  
read_indels_plot(readsetnames, platformlabels)
dev.off()

### Plot the full figure together:
pdf("Figure4Output/Figure4Multiplot.pdf", width=13, height=8)
par(mfrow=c(2,3))
read_mononucqvscore_plot(readsetnames, platformlabels, strtype='mononuc', xlabel=c("Run length"), ymax=25)
read_mononucqvscore_plot(readsetnames, platformlabels, strtype='dinuc', ymax=25, minlength=25, xlabel=c("Run length"), plottitle='Accuracy of dinucleotide runs')
read_mononucqvscore_plot(readsetnames, platformlabels, strtype='trinuc', ymax=25, minlength=25, xlabel=c("Run length"), plottitle='Accuracy trinucleotide runs')
read_qv_plot(readsetnames, platformlabels)
read_substitutions_plot(readsetnames, platformlabels, outputdir="Figure4Output")
read_indels_plot(readsetnames, platformlabels)
dev.off()

# compare Illumina old and new
illuminareadsetnames <- c("illumina_2x250", "illumina_googlepcrplus", "illumina_googlepcrfree", "NIST_onso_2024Q1")
illuminaplatformlabels <- c("2x250 (2016)", "Google PCR Plus", "Google PCR Free", "NIST/PB Onso (2024)")

pdf("IlluminaSubsIndelComparison.pdf", width=20, height=11)
par(mfrow=c(1,2))
read_substitutions_plot(illuminareadsetnames, illuminaplatformlabels, outputdir="Figure4Output", legend=TRUE, legendposx=7)
read_indels_plot(illuminareadsetnames, illuminaplatformlabels, legendposx=7, legendposy=25)
dev.off()

pdf("IlluminaMononucAccuracy.pdf", width=11, height=11)
read_mononucqvscore_plot(illuminareadsetnames, illuminaplatformlabels, errorbars=TRUE, plottitle='Illumina Accuracy of Mononucleotide Runs')
dev.off()

# compare Illumina old and new
shortreadsetnames <- c("illumina_2x250", "illumina_googlepcrplus_ds0.1", "illumina_googlepcrfree_ds0.1", "NIST_onso_2024Q1", "element_ultraq_jun2024")
shortplatformlabels <- c("Illumina 2x250 (2016)", "Illumina PCR Plus (2020)", "Illumina PCR Free (2020)", "NIST/PB Onso (2024)", "Element Ultraq (2024)")
shortreadcolors <- c("#88CCEE", "#AA4499", "#661100", "#44AA99", "#332288")


pdf("ShortReadSubsIndelComparison.pdf", width=23, height=11)
par(mfrow=c(1,2))
read_substitutions_plot(shortreadsetnames, shortplatformlabels, platformcolors=shortreadcolors, outputdir="Figure4Output", legend=TRUE, legendcex=1.2, titlecex=1.3, legendposx=10, legendposy=3700)
read_indels_plot(shortreadsetnames, shortplatformlabels, platformcolors=shortreadcolors, legendcex=1.2, titlecex=1.3, legendposx=10, legendposy=28)
dev.off()

pdf("ShortReadMononucAccuracy.pdf", width=11, height=11)
#read_mononucqvscore_plot(shortreadsetnames, shortplatformlabels, errorbars=TRUE, plottitle='Short Read Accuracy of Mononucleotide Runs', titlecex=1.3, legendcex=1.2)
read_mononucqvscore_plot(shortreadsetnames, shortplatformlabels, errorbars=TRUE, plottitle='', legendcex=1.2)
dev.off()



