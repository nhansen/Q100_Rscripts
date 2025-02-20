setwd("/Users/nhansen/HG002_diploid_benchmark/PaperFigures/CorrectionComparison")

library(colorspace)
library(Hmisc)

################
### FIGURE 4 ###
### Adam's description: ###
### Sequencing coverage and accuracy of various platforms compared against ###'
### the HG002 genome benchmark. [Essential read quality statistics from ###
### Nanopore, PacBio, Illumina, Element, etc. illustrating subtle base ###
### call and coverage biases not evident from typically reported QV scores. ###
### Things like real vs. reported QV values, indel biases, coverage ###
### uniformity, transition/transversion biases ###
################

outdir <- "dorado_vs_herro"
readsetqvfiles <- c("dorado_vs_herro/herro_corrected.errorqvstats.txt", "dorado_vs_herro/dorado_corrected.errorqvstats.txt")
readsetnames <- c("hg002v1.1_dorado_corrected", "hg002v1.1_herro_corrected")
platformlabels <- c("Dorado corrected", "Herro corrected")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#readplatformcolors <- c("#6699CC", "#CC6677", "#01541F", "#332288", "#661100", "#5D5D5D")
readplatformcolors <- c("#6699CC", "#01541F", "#332288", "#661100", "#5D5D5D")
platformlinetype <- c(1, 2, 3, 4, 5, 6)
platformpchvals <- c(0, 1, 2, 5, 6, 8)
filledplatformpchvals <- c(15, 16, 17, 23, 25, 8)

#platformpchvals <- c(15, 17, 19, 3, 2, 0)

title <- ""

qvalues <- function(errorlist, totallist) {
  qscores <- sapply(seq(1, length(errorlist)), function(x) {return(ifelse(totallist[x]==0, NA, as.numeric(-10.0*log10((errorlist[x]+1)/(totallist[x]+1)))))})
  return(qscores)
}

get_qv_counts <- function(file) {
  qvcounts <- read.table(file, header=FALSE, sep="\t")
  names(qvcounts) <- c("QVReported", "SNVErrors", "IndelErrors", "TotalBases")
  
  qvcounts$SNVObsQV <- qvalues(qvcounts$SNVErrors, qvcounts$TotalBases)
  snvaccconfints <- binconf(qvcounts$SNVErrors, qvcounts$TotalBases, return.df=TRUE)
  qvcounts$SNVObsQVHigh <- as.numeric(-10.0*log10(snvaccconfints$Lower+0.0000000001))
  qvcounts$SNVObsQVLow <- as.numeric(-10.0*log10(snvaccconfints$Upper))
  
  qvcounts$IndelObsQV <- qvalues(qvcounts$IndelErrors, qvcounts$TotalBases)
  indelaccconfints <- binconf(qvcounts$IndelErrors, qvcounts$TotalBases, return.df=TRUE)
  qvcounts$IndelObsQVHigh <- as.numeric(-10.0*log10(indelaccconfints$Lower+0.0000000001))
  qvcounts$IndelObsQVLow <- as.numeric(-10.0*log10(indelaccconfints$Upper))
  
  qvcounts$TotalObsQV <- qvalues((qvcounts$SNVErrors+qvcounts$IndelErrors), qvcounts$TotalBases)
  totalaccconfints <- binconf(qvcounts$IndelErrors+qvcounts$SNVErrors, qvcounts$TotalBases, return.df=TRUE)
  qvcounts$TotalObsQVHigh <- as.numeric(-10.0*log10(totalaccconfints$Lower+0.0000000001))
  qvcounts$TotalObsQVLow <- as.numeric(-10.0*log10(totalaccconfints$Upper))
  
  return(qvcounts)
}

read_qv_plot <- function(readsetnames, platformlabels, cexfactor=0.06, plottitle=title, errorbars=FALSE) {
  readsetfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".errorqvstats.txt"), sep="", collapse="")})
  maxqvreported <- 0
  for (i in seq(1, length(readsetfiles))) {
    qvcounts <- get_qv_counts(readsetfiles[i])
    maxreported <- max(qvcounts[qvcounts$TotalBases>0, "QVReported"])
    if (maxreported > maxqvreported) {
      maxqvreported <- maxreported
    }
  }
  
  plot(c(0, maxqvreported), c(0, maxqvreported), type="l", lty=1, main=plottitle, xlab="Reported QV Score", ylab="Observed QV")
  
  for (i in seq(1, length(readsetfiles))) {
    qvcounts <- get_qv_counts(readsetfiles[i])
    qvcounts <- qvcounts[qvcounts$TotalBases>0, ]
    points(qvcounts$QVReported, qvcounts$TotalObsQV, type="b", lty=1, pch=filledplatformpchvals[i], col=readplatformcolors[i], bg=readplatformcolors[i])
    if(errorbars) {
      arrows(x0=qvcounts$QVReported, y0=qvcounts$TotalObsQVLow, x1=qvcounts$QVReported, y1=qvcounts$TotalObsQVHigh, code=3, angle=90, length=0.05, col=readplatformcolors[i])
    }
  }
  legend(5, maxqvreported-5, platformlabels, pch=filledplatformpchvals, lty=platformlinetype, col=readplatformcolors, pt.bg=readplatformcolors)

}

read_qv_density_plot <- function(readsetfiles, platformlabels, cexfactor=0.06, plottitle=title) {
  
  maxqvobserved <- 0
  for (i in seq(1, length(readsetfiles))) {
    qvcounts <- get_qv_counts(readsetfiles[i])
    maxobserved <- max(qvcounts[qvcounts$TotalBases>0, "TotalObsQV"], na.rm=TRUE)
    if (maxobserved > maxqvobserved) {
      maxqvobserved <- maxobserved
    }
  }
  
  firstqvcounts <- get_qv_counts(readsetfiles[1])
  firstqvcounts <- firstqvcounts[order(firstqvcounts$TotalObsQV), ]
  plot(firstqvcounts$TotalObsQV, firstqvcounts$TotalBases/sum(firstqvcounts$TotalBases), type="l", lty=2, main=plottitle, xlab="Benchmark QV Score", ylab="Base density", ylim=c(0,1))
  
  for (i in seq(2, length(readsetfiles))) {
    qvcounts <- get_qv_counts(readsetfiles[i])
    qvcounts <- qvcounts[order(qvcounts$TotalObsQV), ]
    points(qvcounts$TotalObsQV, qvcounts$TotalBases/sum(qvcounts$TotalBases), type="l", lty=1, col=readplatformcolors[i])
  }
  legend("topleft", platformlabels, lty=1, col=readplatformcolors)

}

#pdf("Figure4Output/QVAccuracy.pdf")
#read_qv_plot(readsetnames, platformlabels)
#dev.off()


# Plot mononucleotide accuracy

readmononuchist <- function(filename) {
  mndf <- read.table(filename, header=TRUE, sep="")
  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  mndf <- mndf[mndf$readlength != -1, ]
  mndf$reflength <- as.integer(mndf$reflength)
  mndf$readlength <- as.integer(mndf$readlength)
  mndf$numcorrect <- as.integer(mndf$numcorrect)
  mndf$numerror <- as.integer(mndf$numerror)
  
  return(mndf)
}

plotmononucaccuracy <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Accuracy of mononucleotide runs")
  }
  accarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numcorrect"])/sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})
  #  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  if (addtoplot) {
    points(seq(minlength, maxlength), accarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), ylim=c(0,1), col=color)
  }  
  else {
    plot(seq(minlength, maxlength), accarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), main=plottitle, ylim=c(0,1), col=color)
  }
}

plotmononucerrors <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Error rate in mononucleotide runs")
  }
  errorarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numerror"])/sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})
  #  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  if (addtoplot) {
    points(seq(minlength, maxlength), errorarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Error rate"), ylim=c(0,1), col=color)
  }  
  else {
    plot(seq(minlength, maxlength), errorarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Error rate"), main=plottitle, ylim=c(0,1), col=color)
  }
}

plotmononucqvscores <- function(file, plottitle=NA, titlecex=1.0, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE, ymax=35, errorbars=FALSE){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Accuracy of mononucleotide runs")
  }
  errorarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numerror"])})
  totalarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})
  qvarray <- sapply(seq(minlength, maxlength), function(x) {-10.0*log10(sum(mncounts[mncounts$reflength==x, "numerror"])/sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"]))})
  totalerrorrateconfints <- binconf(errorarray, totalarray, return.df=TRUE)
  qvhigharray <- as.numeric(-10.0*log10(totalerrorrateconfints$Lower+0.0000000001))
  qvlowarray <- as.numeric(-10.0*log10(totalerrorrateconfints$Upper))
  if (addtoplot) {
    points(seq(minlength, maxlength), qvarray, pch=pchval, col=color)
    if (errorbars) {
      arrows(x0=seq(minlength, maxlength), y0=qvlowarray, x1=seq(minlength, maxlength), y1=qvhigharray, code=3, angle=90, length=0.05, col=color)
    }
  }  
  else {
    plot(seq(minlength, maxlength), qvarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Phred QV Score"), main=plottitle, cex.main=titlecex, ylim=c(0,ymax), col=color)
    if (errorbars) {
      arrows(x0=seq(minlength, maxlength), y0=qvlowarray, x1=seq(minlength, maxlength), y1=qvhigharray, code=3, angle=90, length=0.05, col=color)
    }
  }
}

plotmononuccoveragecounts <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE, ymax=35){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Number of traversed mononucleotide runs assessed")
  }
  totalarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})/1000
  if (addtoplot) {
    points(seq(minlength, maxlength), totalarray, pch=pchval, col=color, log='y')
  }  
  else {
    plot(seq(minlength, maxlength), totalarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Thousands of reads"), main=plottitle, col=color, log='y', ylim=c(1,1500))
  }
}
read_mononucacc_plot <- function(readsetnames, platformlabels, maxlength=40, platformcolors=readplatformcolors) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".mononuchist.txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucaccuracy(mnstatsfiles[1], minlength=10, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1])
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucaccuracy(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE)
    }
  }
  legend("bottomleft", platformlabels, col=platformcolors, pch=platformpchvals)
}

read_mononucerror_plot <- function(readsetnames, platformlabels, maxlength=40, platformcolors=readplatformcolors) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".mononuchist.txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucerrors(mnstatsfiles[1], minlength=10, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1])
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucerrors(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE)
    }
  }
  legend(30, 0.4, platformlabels, col=platformcolors, pch=platformpchvals)
}

read_mononucqvscore_plot <- function(readsetnames, platformlabels, maxlength=40, ymax=45, platformcolors=readplatformcolors, errorbars=FALSE, plottitle=NA, titlecex=1.0, legendcex=1.0) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".mononuchist.txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucqvscores(mnstatsfiles[1], minlength=10, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1], ymax=ymax, errorbars=errorbars, plottitle=plottitle, titlecex=titlecex)
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucqvscores(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE, ymax=ymax, errorbars=errorbars)
    }
  }
  legend(32, 42, platformlabels, col=platformcolors, pch=platformpchvals, cex=legendcex)
}

read_mononuccoverage_plot <- function(mnstatsfiles, platformlabels, maxlength=40, ymax=45, platformcolors=readplatformcolors) {
  
  if (length(mnstatsfiles)>0) {
    plotmononuccoveragecounts(mnstatsfiles[1], minlength=10, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1], ymax=ymax)
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononuccoveragecounts(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE, ymax=ymax)
    }
  }
  legend(32, 500, platformlabels, col=platformcolors, pch=platformpchvals)
}
### Make mononucleotide accuracy plot ###


pdf("dorado_vs_herro/ReadMononucAccuracy.pdf")  
read_mononucacc_plot(readsetnames, platformlabels)
dev.off()

pdf("dorado_vs_herro/ReadMononucErrorRate.pdf")  
read_mononucerror_plot(readsetnames, platformlabels)
dev.off()

pdf("dorado_vs_herro/ReadMononucQVScores.pdf")  
read_mononucqvscore_plot(readsetnames, platformlabels)
dev.off()
# Plot substitution rate-by-type histogram

typeorder <- c("A_C", "A_G", "A_T", "T_C", "T_G", "T_A", "G_A", "G_T", "G_C", "C_A", "C_T", "C_G")
titv <- c("tv", "ti", "tv", "ti", "tv", "tv", "ti", "tv", "tv", "tv", "ti", "tv" )

read_substitutions_plot <- function(readsetnames, platformlabels, platformcolors=readplatformcolors, outputdir, xlabval="Read platform", ylabval="Substitutions per mb", titleval="Read substitution error rates", titlecex=1.0, ymax=NA, legend=TRUE, legendposx=NA, legendposy=NA, legendcex=1) {
  subsfilenames <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".singlenucerrorstats.txt"), sep="", collapse="")})
  
  firsthist <- read.table(subsfilenames[1], sep="\t")
  names(firsthist) <- c("errortype", "errorcount", "errorspermbaligned")
  typeorderindex <- sapply(firsthist$errortype, function(x) {which(typeorder==x)})
  firsttis <- sum(firsthist[titv[typeorderindex]=="ti", "errorspermbaligned"])
  firsttvs <- sum(firsthist[titv[typeorderindex]=="tv", "errorspermbaligned"])
  tivals <- c(firsttis)
  tvvals <- c(firsttvs)
  platformlabelswithtitv <- sapply(platformlabels, function(x) {label = paste(c("Transitions       Transversions\n", x), sep="", collapse=""); return(label)})
  
  for (i in seq(2, length(subsfilenames))) {
    subshist <- read.table(subsfilenames[i], sep="\t")
    names(subshist) <- c("errortype", "errorcount", "errorspermbaligned")
    typeorderindex <- sapply(subshist$errortype, function(x) {which(typeorder==x)})
    
    platformtis <- sum(subshist[titv[typeorderindex]=="ti", "errorspermbaligned"])
    platformtvs <- sum(subshist[titv[typeorderindex]=="tv", "errorspermbaligned"])
    
    tivals <- append(tivals, platformtis)
    tvvals <- append(tvvals, platformtvs)
  }
  
  barcolors <- sapply(platformcolors, function(x) {c(darken(x), lighten(x, 0.3))})
  
  if (is.na(ymax)) {
    out <- barplot(rbind(tivals, tvvals), names.arg=platformlabelswithtitv, beside=TRUE, col=barcolors, main=titleval, cex.main=titlecex, xlab=xlabval, ylab=ylabval)
  }
  else {
    out <- barplot(rbind(tivals, tvvals), names.arg=platformlabelswithtitv, beside=TRUE, col=barcolors, main=titleval, cex.main=titlecex, xlab=xlabval, ylab=ylabval, ylim=c(0,ymax))
  }
  if (is.na(legendposx)) {
    legendposx <- out[2*length(tivals)-4]
  }
  if (is.na(legendposy)) {
    legendposy <- max(rbind(tivals, tvvals))-200
  }
  if (legend) {
    legend(legendposx, legendposy, platformlabels, col=platformcolors, pch=15, cex=legendcex)
  }
  #if (legend) {
    #platformtitvlabels <- as.vector(sapply(platformlabels, function(p) {sapply(c('Ti', 'Tv'), function(titv) {paste(p, titv, sep=" ")})}))
    #legend(legendposx, legendposy, platformtitvlabels, col=barcolors, pch=15, cex=legendcex)
  #}
  formattedrates <- as.integer((rbind(tivals, tvvals)+0.5)*10)/10
  formattedrates <- ifelse(formattedrates>10, as.integer(formattedrates), formattedrates)
  text(out, rbind(tivals, tvvals), formattedrates, pos=3, xpd=NA, cex=0.85)
  
}

pdf("dorado_vs_herro/ReadSubstitutionRates.pdf")
read_substitutions_plot(readsetnames, platformlabels, outputdir="dorado_vs_herro", legend=TRUE, legendposx=4, legendposy=170)
dev.off()


### Indel error plot ###

totalindelrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist$indelspermbaligned, na.rm=TRUE))  
}

totalinsertionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist[indellengthhist$indellength>0, "indelspermbaligned"], na.rm=TRUE))  
}

totaldeletionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist[indellengthhist$indellength<0, "indelspermbaligned"], na.rm=TRUE))  
}

read_indels_plot <- function(indelreadsetnames, platformlabels, platformcolors=readplatformcolors, xlabval="Read platform", ylabval="Indels per mb", titleval="Indel errors in reads", titlecex=1.0, ymax=NA, legendcex=1.0, legendposx=NA, legendposy=NA) {
  indellengthfiles <- sapply(indelreadsetnames, function(x) {paste(c(outdir, "/", x, ".indelerrorstats.txt"), sep="", collapse="")})
  totalinsertionrates <- sapply(indellengthfiles, totalinsertionrate)
  totaldeletionrates <- sapply(indellengthfiles, totaldeletionrate)
  
  platformlabelswithinsdel <- sapply(platformlabels, function(x) {label = paste(c("Insertions       Deletions\n", x), sep="", collapse=""); return(label)})
  
  barcolors <- sapply(platformcolors, function(x) {c(darken(x), lighten(x, 0.3))})
  out <- barplot(rbind(totalinsertionrates, totaldeletionrates), beside=TRUE, names.arg=platformlabelswithinsdel, col=barcolors, main=titleval, cex.main=titlecex, xlab=xlabval, ylab=ylabval)
  if (is.na(legendposx)) {
    legendposx <- out[2*length(totalinsertionrates)-2]
  }
  if (is.na(legendposy)) {
    legendposy <- max(rbind(totalinsertionrates, totaldeletionrates))-100
  }
  legend(legendposx, legendposy, platformlabels, col=platformcolors, pch=15, cex=legendcex)
  formattedrates <- as.integer((rbind(totalinsertionrates, totaldeletionrates)+0.5)*10)/10
  formattedrates <- ifelse(formattedrates>100, as.integer(formattedrates), formattedrates)
  text(out, rbind(totalinsertionrates, totaldeletionrates), formattedrates, pos=3, xpd=NA, cex=0.85)
}

pdf("dorado_vs_herro/ReadIndelRates.pdf")  
read_indels_plot(readsetnames, platformlabels, legendposx=4, legendposy=500)
dev.off()

### Plot the full figure together:
pdf("Figure4Multiplot.pdf", width=11, height=11)
par(mfrow=c(2,2))
read_qv_plot(readsetqvfiles, platformlabels)
read_mononucqvscore_plot(readsetnames, platformlabels, errorbars=TRUE)
read_substitutions_plot(readsetnames, platformlabels, outputdir="Figure4Output")
read_indels_plot(readsetnames, platformlabels)
dev.off()

# compare Illumina old and new
illuminareadsetnames <- c("illumina_2x250", "illumina_googlepcrplus_ds0.1", "illumina_googlepcrfree_ds0.1", "NIST_onso_2024Q1")
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



