setwd("/Users/nhansen/HG002_diploid_benchmark/plots/ngaplots")

args = commandArgs(trailingOnly=TRUE)
assemblycolors <- c("#44AA99", "#332288", "#882255", "#888888", "#DDCC77", "#661100", "#6699CC")
methodcolors <- c("#DDCC77", "#661100", "#6699CC")

assemblyname <- ifelse(!( is.na(args[1])), args[1], "assembly")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c("NGAx ", assemblyname, " vs ", genomename), sep="", collapse=""))

assemblylabels <- c("Q28/HiFi Hap1", "Q28/HiFi Hap2", "Q28/bothraw/def Hap1", "Q28/bothraw/def Hap2", "Q28/bothraw/hifiasm Hap1", "Q28/bothraw/hifiasm Hap2")
assemblysizefiles <- c("q28_hifi_hap1.alignclusterlengths.txt", "q28_hifi_hap2.alignclusterlengths.txt", "q28_both_raw_hap1.alignclusterlengths.txt", "q28_both_raw_hap2.alignclusterlengths.txt", "q28_bothraw_hifiasm_hap1.alignclusterlengths.txt", "q28_bothraw_hifiasm_hap2.alignclusterlengths.txt")
#assemblysizefiles <- sapply(assemblysizefiles, function(x) {file=paste(c(outputdir, "/", x), sep="", collapse=""); return(file)})

readlengths <- function(sizefile) {
  clusterlengths <- read.table(sizefile, sep="\t", header=FALSE)
  names(clusterlengths) <- c("perc", "clusterlength", "totallength", "chromosome")
  clusterlengths$clusterlength <- clusterlengths$clusterlength/1000000
  
  return(clusterlengths)
}

plotclusterlengths <- function(clusterlengths, color="red", title="NGAx", dashed=FALSE, ltyval=NA, cexval=1.0) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  if (is.na(ltyval)) {
    ltyval <- ifelse(dashed, 2, 1)
  }
  plot(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,100), ylim=c(0,250), xlab="Percent of Haploid Genome", ylab="Aligned length", main=title, cex=cexval)
}

addclusterlengths <- function(clusterlengths, color="blue", dashed=FALSE, ltyval=NA) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  if (is.na(ltyval)) {
    ltyval <- ifelse(dashed, 2, 1)
  }
  lines(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,100))
}


assembly_ngax_plot <- function(clusterfiles, contigfiles=c(), scaffoldfiles=c(), assemblylabels=c(), ideal=FALSE, haplotype="MAT", plottitle="", cexval=1.0) {
  
  firstclusters <- readlengths(clusterfiles[1]) 
  plotclusterlengths(firstclusters, col=assemblycolors[1], title=plottitle, lty=1, cexval=cexval)
  if (length(contigfiles) > 1) {
    firstcontigs <- readlengths(contigfiles[1]) 
    addclusterlengths(firstcontigs, col=assemblycolors[1])    
  }
  if (length(scaffoldfiles) > 1) {
    firstscaffolds <- readlengths(scaffoldfiles[1]) 
    addclusterlengths(firstscaffolds, col=assemblycolors[1])    
  }
  
  if (length(clusterfiles) > 1) {
    for (i in seq(2, length(clusterfiles))) {
      #addclusterlengths(readlengths(clusterfiles[i]), col=assemblycolors[i], lty=i)
      addclusterlengths(readlengths(clusterfiles[i]), col=assemblycolors[i], lty=1)
    }
  }
  if (ideal) {
    ideallengths <- readlengths("Figure3Output/v1.1.mat.alignclusterlengths.txt")
    addclusterlengths(ideallengths, col="black")
  }
  #legend("topright", assemblylabels, col=assemblycolors, lty=seq(1, length(clusterfiles)))
  legend("topright", assemblylabels, col=assemblycolors, lty=rep(1, length(clusterfiles)))
  
}

# Make NGAx plot:

pdf("NGAxPlot.pdf")
assembly_ngax_plot(assemblysizefiles, assemblylabels=assemblylabels, ideal=FALSE, plottitle="NGAx for different assemblies")
dev.off()
