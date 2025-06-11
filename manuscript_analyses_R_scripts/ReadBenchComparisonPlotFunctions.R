library(Hmisc)

# function get_qv_counts reads a file of QV counts produced by GQC's readbench and returns a dataframe
# of binned reported QV scores with columns:
# QVReported SNVErrors IndelErrors TotalBases SNVObsQV SNVObsQVLow SNVObsQVHigh 
# IndelObsQV IndelObsQVLow IndelObsQVHigh TotalObsQV TotalObsQVLow TotalObsQVHigh
# where the observed QV values are the binomial expectation value and its Bayesian
# confidence interval (calculated by the Hmisc library's "binconf" function)

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

qvalues <- function(errorlist, totallist) {
  qscores <- sapply(seq(1, length(errorlist)), function(x) {return(ifelse(totallist[x]==0, NA, as.numeric(-10.0*log10((errorlist[x]+1)/(totallist[x]+1)))))})
  return(qscores)
}

# Plot QV statistics: so long as a qvstats file exists at the location outdir/readsetname.errorqvstats.txt, 
# the read_qv_plot function plots an x/y dot plot of observed vs. reported QV score

read_qv_plot <- function(readsetnames, platformlabels, cexfactor=0.06, plottitle="Read QV score accuracy", errorbars=FALSE, legend=TRUE, legendx=NA, legendy=NA) {
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
  if (legend) {
    if (is.na(legendx)) {
      legendx = 5
    }
    if (is.na(legendy)) {
      legendy = maxqvreported-5
    }  
    legend(legendx, legendy, platformlabels, pch=filledplatformpchvals, lty=platformlinetype, col=readplatformcolors, pt.bg=readplatformcolors)
  }
  
}

read_qv_density_plot <- function(readsetnames, platformlabels, cexfactor=0.06, plottitle="Fraction of read bases for observed QVs") {
  readsetfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".errorqvstats.txt"), sep="", collapse="")})
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
  plot(firstqvcounts$TotalObsQV, firstqvcounts$TotalBases/sum(firstqvcounts$TotalBases), type="l", lty=1, col=readplatformcolors[1], main=plottitle, xlab="Benchmark QV Score", ylab="Base density", ylim=c(0,1))
  
  for (i in seq(2, length(readsetfiles))) {
    qvcounts <- get_qv_counts(readsetfiles[i])
    qvcounts <- qvcounts[order(qvcounts$TotalObsQV), ]
    points(qvcounts$TotalObsQV, qvcounts$TotalBases/sum(qvcounts$TotalBases), type="l", lty=1, col=readplatformcolors[i])
  }
  legend("topleft", platformlabels, lty=1, col=readplatformcolors)
}

# Plot substitution rate-by-type histogram: so long as a snvstats file exists at the location 
# outdir/readsetname.singlenucerrorstats.txt, with three tab-delimited columns containing the
# "typeorder" value (see variable below), the total number of errors, and the number of errors
# per total Mb of read bases (not just per total number of that particular reference base),
# this routine will produce a histogram of substitution rates divided into transitions and
# transversions

typeorder <- c("A_C", "A_G", "A_T", "T_C", "T_G", "T_A", "G_A", "G_T", "G_C", "C_A", "C_T", "C_G")
titv <-      c("tv",  "ti",  "tv",  "ti",  "tv",  "tv",  "ti",  "tv",  "tv",  "tv",  "ti",  "tv" )

read_substitutions_plot <- function(readsetnames, platformlabels, platformcolors=readplatformcolors, outputdir, xlabval="Read platform", ylabval="Substitutions per mb", titleval="Read substitution error rates", titlecex=1.0, ymax=NA, legend=TRUE, legendposx=NA, legendposy=NA, legendcex=1) {
  subsfilenames <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".singlenucerrorstats.txt"), sep="", collapse="")})
  
  firsthist <- read.table(subsfilenames[1], sep="\t")
  names(firsthist) <- c("errortype", "errorcount", "errorspermbaligned")
  typeorderindex <- sapply(firsthist$errortype, function(x) {which(typeorder==x)})
  firsttis <- sum(firsthist[titv[typeorderindex]=="ti", "errorspermbaligned"])
  firsttvs <- sum(firsthist[titv[typeorderindex]=="tv", "errorspermbaligned"])
  tivals <- c(firsttis)
  tvvals <- c(firsttvs)
  platformlabelswithtitv <- sapply(platformlabels, function(x) {label = paste(c("Ti    Tv\n", x), sep="", collapse=""); return(label)})
  
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

### Indel error plot: so long as an indel error stats file exists at the location 
# outdir/readsetname.indelerrorstats.txt, with three tab-delimited columns containing the
# difference between number of bases in the reads and number of bases in the benchmark
# (positive values are insertions), the total number of errors, and the number of errors
# per total Mb of read bases (not just per total number of that particular reference base),
# the read_indels_plot function will produce a histogram of indel rates divided into 
# insertions and deletions
###

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
  
  platformlabelswithinsdel <- sapply(platformlabels, function(x) {label = paste(c("Ins   Del\n", x), sep="", collapse=""); return(label)})
  
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

# Plot mononucleotide accuracy

readmononuchist <- function(filename) {
  mndf <- read.table(filename, header=TRUE, sep="\t")
  names(mndf) <- c("reflength", "numcorrect", "numlengtherrors", "numcomplexerrors", "numflankerrors", "numhetalleles")
  
  return(mndf)
}

plotmononucaccuracy <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE) {
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
  validreflengths <- mncounts[which(mncounts$reflength>=minlength & mncounts$reflength<=maxlength), "reflength"]
  accarray <- sapply(validreflengths, function(x) { return(mncounts[mncounts$reflength==x, "numcorrect"]/(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"]))})
  if (addtoplot) {
    points(validreflengths, accarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), ylim=c(0,1), col=color)
  }  
  else {
    plot(validreflengths, accarray, pch=pchval, xlab=c("Mononucleotide run length"), xlim=c(minlength, maxlength), ylab=c("Accuracy"), main=plottitle, ylim=c(0,1), col=color)
  }
}

plotmononucerrors <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE) {
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
  validreflengths <- mncounts[which(mncounts$reflength>=minlength & mncounts$reflength<=maxlength), "reflength"]
  errorarray <- sapply(validreflengths, function(x) {return(1.0-mncounts[mncounts$reflength==x, "numcorrect"]/(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"]))})
  #  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  if (addtoplot) {
    points(validreflengths, errorarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Error rate"), ylim=c(0,1), col=color)
  }  
  else {
    plot(validreflengths, errorarray, pch=pchval, xlab=c("Mononucleotide run length"), xlim=c(minlength, maxlength), ylab=c("Error rate"), main=plottitle, ylim=c(0,1), col=color)
  }
}

plotmononucqvscores <- function(file, plottitle=NA, titlecex=1.0, minlength=NA, maxlength=NA, pchval=16, color="black", xlabel=c("Mononucleotide run length"), addtoplot=FALSE, ymax=35, errorbars=FALSE){
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
  validreflengths <- mncounts[which(mncounts$reflength>=minlength & mncounts$reflength<=maxlength), "reflength"]
  errorarray <- sapply(validreflengths, function(x) {return(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"])})
  totalarray <- sapply(validreflengths, function(x) {return(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"])})
  #qvarray <- sapply(validreflengths, function(x) {-10.0*log10(mncounts[mncounts$reflength==x, "numerror"]/(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"]))})
  qvarray <- -10.0*log10(errorarray/totalarray)
  totalerrorrateconfints <- binconf(errorarray, totalarray, return.df=TRUE)
  qvhigharray <- as.numeric(-10.0*log10(totalerrorrateconfints$Lower+0.0000000001))
  qvlowarray <- as.numeric(-10.0*log10(totalerrorrateconfints$Upper))
  if (addtoplot) {
    points(validreflengths, qvarray, pch=pchval, col=color)
    if (errorbars) {
      arrows(x0=validreflengths, y0=qvlowarray, x1=validreflengths, y1=qvhigharray, code=3, angle=90, length=0.05, col=color)
    }
  }  
  else {
    plot(validreflengths, qvarray, pch=pchval, xlab=xlabel, ylab=c("Phred QV Score"), main=plottitle, cex.main=titlecex, ylim=c(0,ymax), col=color)
    if (errorbars) {
      arrows(x0=validreflengths, y0=qvlowarray, x1=validreflengths, y1=qvhigharray, code=3, angle=90, length=0.05, col=color)
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
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".strlengthaccuracy.mononuc.txt"), sep="", collapse="")})
  
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

read_mononucerror_plot <- function(readsetnames, platformlabels, minlength=10, maxlength=40, strtype='mononuc', platformcolors=readplatformcolors) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".strlengthaccuracy.", strtype, ".txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucerrors(mnstatsfiles[1], minlength=minlength, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1])
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucerrors(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE)
    }
  }
  legend(30, 0.4, platformlabels, col=platformcolors, pch=platformpchvals)
}

read_mononucqvscore_plot <- function(readsetnames, platformlabels, minlength=NA, maxlength=40, xlabel=c("Mononucleotide run length"), ymax=45, strtype='mononuc', platformcolors=readplatformcolors, errorbars=FALSE, plottitle=NA, titlecex=1.0, legendcex=1.0) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".strlengthaccuracy.", strtype, ".txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucqvscores(mnstatsfiles[1], minlength=10, maxlength=maxlength, xlabel=xlabel, pchval=platformpchvals[1], color=platformcolors[1], ymax=ymax, errorbars=errorbars, plottitle=plottitle, titlecex=titlecex)
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucqvscores(mnstatsfiles[i], maxlength=maxlength, xlabel=xlabel, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE, ymax=ymax, errorbars=errorbars)
    }
  }
  legend(32, 42, platformlabels, col=platformcolors, pch=platformpchvals, cex=legendcex)
}

read_mononuccoverage_plot <- function(mnstatsfiles, platformlabels, maxlength=40, ymax=45, xlabel=c("Mononucleotide run length"), platformcolors=readplatformcolors) {
  
  if (length(mnstatsfiles)>0) {
    plotmononuccoveragecounts(mnstatsfiles[1], minlength=10, maxlength=maxlength, xlabel=xlabel, pchval=platformpchvals[1], color=platformcolors[1], ymax=ymax)
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononuccoveragecounts(mnstatsfiles[i], maxlength=maxlength, xlabel=xlabel, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE, ymax=ymax)
    }
  }
  legend(32, 500, platformlabels, col=platformcolors, pch=platformpchvals)
}
