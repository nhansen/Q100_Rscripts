setwd("/Users/nhansen/HG002_diploid_benchmark/PaperFigures")

################
### FIGURE 1 ###
### Adam's description: ###
### Overview of the complete diploid HG002 genome. ###
### [Diploid ideogram of all HG002 chromosomes, with maternal/paternal homologs ###
### adjacent to one another and each chromosome labeled with things like ###
### satellite annotation, het SNV density, het SVs, suspicious regions, gaps, etc.] ###
################

library(stringr)
library(karyoploteR)
library(ggplot2)
library(grid)
library(gridExtra)
library("cowplot")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

strandcolors <- c("#332288", "#88CCEE")

censat_color_palette <- c("#662F90",
                          "#FA99FF",
                          "#AC33C7",
                          "#00CCCC",
                          "#990000",
                          "#FF6600",
                          "#FF9200",
                          "#FFCC99",
                          "#00DE60",
                          "#1B998B",
                          "#0080FA",
                          "#335189",
                          "#E0E0E0",
                          "#000000",
                          "#FFCC00",
                          "#CC0000",
                          "#78A1BB",
                          "#FF6600",
                          "#FFCC99",
                          "#990000")
censat_labels <- c("rDNA",
  "bsat",
  "gsat",
  "censat",
  "activeHOR",
  "inactiveHOR",
  "dHOR",
  "mon",
  "HSat1A",
  "HSat1B",
  "HSat2",
  "HSat3",
  "ct",
  "GAP",
  "subTerm",
  "mixedAlpha",
  "HSat2_3", 
  "hor(SF3-5)",
  "mon/hor(SF4)",
  #"active_hor(SF1, SF01, SF4)"
  "active_hor")

genomefile <- "KaryoploteRFiles/v1.1.karyotype.txt"
matgenomefile <- "KaryoploteRFiles/v1.1.karyotype.mat.txt"
patgenomefile <- "KaryoploteRFiles/v1.1.karyotype.pat.txt"
benchgenome <- toGRanges(genomefile)
matbenchgenome <- toGRanges(matgenomefile)
patbenchgenome <- toGRanges(patgenomefile)

#censatfile <- "KaryoploteRFiles/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed"
#censats <- toGRanges(censatfile)

matcensatfile <- "KaryoploteRFiles/hg002v1.1.cenSatv2.0.MAT.bed"
matcensats <- toGRanges(matcensatfile)

patcensatfile <- "KaryoploteRFiles/hg002v1.1.cenSatv2.0.PAT.bed"
patcensats <- toGRanges(patcensatfile)

allcensats <- c(matcensats, patcensats)

hetfile <- "KaryoploteRFiles/v1.1.hets1mb.1to1alignable.windows.withhetcounts.norm50.heatmap.bed"
heterozygosity <- toGRanges(hetfile)

non1to1hetfile <- "KaryoploteRFiles/v1.1.hets.non1to1.200k.hapaligned.windows.withheatmap.bed"
non1to1heterozygosity <- toGRanges(non1to1hetfile)

matnon1to1hetfile <- "KaryoploteRFiles/v1.1.hets.non1to1alignable.windows.withhetcounts.heatmap.mat.bed"
patnon1to1hetfile <- "KaryoploteRFiles/v1.1.hets.non1to1alignable.windows.withhetcounts.heatmap.pat.bed"
matnon1to1heterozygosity <- toGRanges(matnon1to1hetfile)
patnon1to1heterozygosity <- toGRanges(patnon1to1hetfile)

excludedfile <- "KaryoploteRFiles/v1.1.excluded_regions.bed"
excluded <- toGRanges(excludedfile)

matexcludedfile <- "KaryoploteRFiles/v1.1.excluded_regions.mat.bed"
matexcluded <- toGRanges(matexcludedfile)

patexcludedfile <- "KaryoploteRFiles/v1.1.excluded_regions.pat.bed"
patexcluded <- toGRanges(patexcludedfile)

plot_homologue_pair <- function(genome, chrom, het=heterozygosity, matcentro=matcensats, patcentro=patcensats, centrolow=0.2, centrohigh=0.4, matloqual=matexcluded, patloqual=patexcluded, loquallow=0.1, loqualhigh=0.3, invert=FALSE) {
  regexstring <- paste(c("^", chrom, "_"), sep="", collapse="")
  allchroms <- seqlevels(genome)
  if (chrom=="chrX") {
    chroms <- c("chrX_MATERNAL", "chrY_PATERNAL")
  }
  else {
    chroms <- allchroms[grep(regexstring, allchroms)]
  }

  pp <- getDefaultPlotParams(plot.type=2)
  pp$ideogramheight <- 400
  pp$data1height <- 1200
  pp$data2height <- 1200
  pp$topmargin <- 650

  kp <- plotKaryotype(genome=genome, plot.type=2, chromosomes=chroms, lwd=0.0, ideogram.plotter=NULL, labels.plotter=NULL, cex=0.4, plot.params=pp)
  kpmidfunction <- kp$ideogram.mid
  offset <- 0
  offsetfunction <- function(chr) {
    return(ifelse(grepl('MATERNAL', chr), kpmidfunction(chr) + offset, kpmidfunction(chr) - offset))
  }
  kp$ideogram.mid <- offsetfunction
  kpAddCytobands(kp, lwd=0.0, labels.plotter=NULL, cex=0.4)

  if (invert) {
    kpRect(kp, data=allcensats, col=allcensats$itemRgb, data.panel="ideogram", lwd=0.0, y0=rep(0, length(allcensats)), y1=rep(1, length(allcensats)))
    
    kpRect(kp, data=het, col=het$itemRgb, data.panel=1, lwd=0.0, y0=rep(centrolow, length(het)), y1=rep(centrohigh, length(het)))
  }
  else {
    kpRect(kp, data=het, col=het$itemRgb, data.panel="ideogram", lwd=0.0, y0=rep(0, length(het)), y1=rep(1, length(het)))
    
    kpRect(kp, data=matcentro, col=matcentro$itemRgb, data.panel=1, lwd=0.0, y0=rep(centrolow, length(matcentro)), y1=rep(centrohigh, length(matcentro)))
    kpRect(kp, data=patcentro, col=patcentro$itemRgb, data.panel=2, lwd=0.0, y0=rep(centrolow, length(patcentro)), y1=rep(centrohigh, length(patcentro)))
  }
  mtext(chrom, side = 2, outer = FALSE)
  # not including low qual regions for now:
  #kpRect(kp, data=patloqual, col="black", data.panel=2, lwd=1.0, y0=rep(loquallow, length(patloqual)), y1=rep(loqualhigh, length(patloqual)))
  #kpRect(kp, data=matloqual, col="black", data.panel=1, lwd=1.0, y0=rep(loquallow, length(matloqual)), y1=rep(loqualhigh, length(matloqual)))
}

plot_all_pairs <- function(genome, het=heterozygosity, mathet=matnon1to1heterozygosity, pathet=patnon1to1heterozygosity, matcentro=matcensats, patcentro=patcensats, centrolow=0.2, centrohigh=0.4, matloqual=matexcluded, patloqual=patexcluded, loquallow=0.1, loqualhigh=0.4, invert=FALSE, labelsat=NA, labelwidths=NA, chromlabels=NA, plottitle="hg002v1.1") {
  
  allchroms <- rev(seqlevels(genome))
  pp <- getDefaultPlotParams(plot.type=2)
  pp$ideogramheight <- 800
  pp$data1height <- 1400
  pp$data2height <- 1400
  pp$topmargin <- 650
  
  kp <- plotKaryotype(genome=genome, chromosomes=allchroms, plot.type=2, lwd=0.0, ideogram.plotter=NULL, labels.plotter=NULL, cex=0.4, plot.params=pp)

  if (invert) {
    #kpAddCytobands(kp, lwd=0.0, labels.plotter=NULL, cex=0.4)
    kpRect(kp, data=matbenchgenome, col="gray", data.panel=1, lwd=0.0, y0=rep(centrolow, length(matbenchgenome)), y1=rep(centrohigh, length(matbenchgenome)))
    kpRect(kp, data=patbenchgenome, col="gray", data.panel=2, lwd=0.0, y0=rep(centrolow, length(patbenchgenome)), y1=rep(centrohigh, length(patbenchgenome)))
    kpRect(kp, data=allcensats, col=allcensats$itemRgb, data.panel="ideogram", lwd=0.0, y0=rep(0, length(allcensats)), y1=rep(1, length(allcensats)))
    kpRect(kp, data=mathet, col=mathet$itemRgb, data.panel=1, lwd=0.0, y0=rep(centrolow, length(matcentro)), y1=rep(centrohigh, length(matcentro)))
    kpRect(kp, data=pathet, col=pathet$itemRgb, data.panel=2, lwd=0.0, y0=rep(centrolow, length(patcentro)), y1=rep(centrohigh, length(patcentro)))
    #mtext("hg002v1.1", side = 2, outer = TRUE)
    #mtext(paste(chroms_to_plot, sep=" ", collapse="   "), side=2, outer=FALSE, cex=0.4)
  }
  else {
    kpAddCytobands(kp, lwd=0.0, labels.plotter=NULL, cex=0.4)
    #kpRect(kp, data=genome, col="white", data.panel="ideogram", lwd=0.0, y0=rep(0, length(het)), y1=rep(1, length(het)))
    kpRect(kp, data=het, col=het$itemRgb, data.panel="ideogram", lwd=0.0, y0=rep(0, length(het)), y1=rep(1, length(het)))
    
    kpRect(kp, data=matcentro, col=matcentro$itemRgb, data.panel=1, lwd=0.0, y0=rep(centrolow, length(matcentro)), y1=rep(centrohigh, length(matcentro)))
    kpRect(kp, data=patcentro, col=patcentro$itemRgb, data.panel=2, lwd=0.0, y0=rep(centrolow, length(patcentro)), y1=rep(centrohigh, length(patcentro)))
    #mtext("hg002v1.1", side = 2, outer = TRUE)
    #mtext(paste(chroms_to_plot, sep=" ", collapse="      "), side=2, line=1, outer=FALSE, cex=0.4)
  }
  if (!is.na(labelsat) && !is.na(labelwidths)) {
     labelpositions <- seq(0, length(chromlabels)-1)*labelwidths + labelsat
  }
  if (length(chromlabels)>1 || !is.na(chromlabels)) {
    mtext(chromlabels, side=2, at=labelpositions, line=1, outer=FALSE, cex=0.4)
  }
  mtext(plottitle, side = 2, outer = TRUE, line=-1)
}

plot_karyoplot_legends <- function() {
  #plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  #legend("topleft", legend =censat_labels, pch=16, pt.cex=3, cex=1.5, bty='n',
         #col = censat_color_palette)
  #mtext("Centromere Sequences", at=0.2, cex=2)

  dataframe <- data.frame(X=rep(0, length(censat_labels)), Y=rep(1, length(censat_labels)), CentromereSequence=factor(censat_labels, levels = censat_labels), itemRgb=censat_color_palette)
  gplot <- ggplot(dataframe, aes(X, Y, fill=CentromereSequence)) + geom_point(size = 7, shape=22) + scale_fill_manual(values=censat_color_palette) + scale_color_manual(values=censat_color_palette) + guides(fill=guide_legend(ncol=4))
  legend <- get_legend(gplot)                     
  
  # Create new plot window 
  grid.newpage()                               
  
  # Draw Only legend  
  grid.draw(legend)  
}

plot_het_heatmap_legend <- function() {
  
  legend(x="right", legend=c("min", "med", "max"),fill=heat.colors(3))
}

chroms_to_plot <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")
ideogramchromlabels <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX/Y")

for (i in seq(length(chroms_to_plot))) {
  chrom <- chroms_to_plot[[i]]
  plotpdffile <- paste("Figure1Output/", chrom, ".figure1karyoplot.pdf", sep="", collapse="")
  pdf(plotpdffile)
  plot_homologue_pair(benchgenome, chrom)
  dev.off()
}

plot_all_pairs(benchgenome)
plot_all_pairs(benchgenome, het=non1to1heterozygosity)

# Can we make a neat plot of the chromosome 8 inverted region?

all_aligns <- "KaryoploteRFiles/all_interhapaligns.paf"
all_align_df <- read.table(all_aligns, sep="\t", header=FALSE)
names(all_align_df) <- c("matchrom", "matchromlength", "matstart", "matend", "strand", "patchrom", "patchromlength", "patstart", "patend", "alignmatches", "totalalignlength", "alignscore", "tag1", "tag2", "tag3", "tag4")
longall_align_df <- all_align_df[all_align_df$totalalignlength>=100000, ]
plot_chromosome_alignments <- function(genome, chromname) {
  regexstring <- paste(c("^", chromname, "_"), sep="", collapse="")
  allchroms <- seqlevels(genome)
  if (chromname=="chrX") {
    chroms <- c("chrX_MATERNAL", "chrY_PATERNAL")
  }
  else {
    chroms <- allchroms[grep(regexstring, allchroms)]
  }
  
  pp <- getDefaultPlotParams(plot.type=1)
  pp$ideogramheight <- 0
  pp$data1height <- 0
  pp$data2height <- 0
  pp$data1inmargin <- 0
  
  chrom_align_df <- longall_align_df[grep(regexstring, longall_align_df$matchrom), ] 

  kp <- plotKaryotype(genome=genome, chromosomes=chroms, plot.type=1, cex=0.4, plot.params=pp)
  kpAddBaseNumbers(kp, tick.dist = 20000000, cex=0.3, minor.ticks=FALSE, tick.len=1)
  inv1start <-toGRanges(chrom_align_df[,c("matchrom", "matstart", "matend")])
  inv1end <- toGRanges(chrom_align_df[,c("patchrom", "patstart", "patend")])
  strand(inv1end) <- chrom_align_df$strand
  forwardcolor <- paste(strandcolors[2], 'CC', sep="", collapse="")
  reversecolor <- paste(strandcolors[1], 'CC', sep="", collapse="")
  aligncolors <- ifelse(chrom_align_df$strand=="+", forwardcolor, reversecolor)
  kpPlotLinks(kp, data=inv1start, data2=inv1end, col=aligncolors, r0=0, r1=0)
}

plot_chromosome_alignments(benchgenome, "chr8")

