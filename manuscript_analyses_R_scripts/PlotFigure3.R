setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/PaperFigures/Figures")

library(colorspace)
library(Hmisc)
library(plotrix)

source("/Users/nhansen/OneDrive/HG002_diploid_benchmark/Q100_Rscripts/manuscript_analyses_R_scripts/AssemblyBenchComparisonPlotFunctions.R")

################
### FIGURE 3 ###
### Adam's description: ###
### Historical assemblies evaluated against the HG002 genome benchmark. ###
### [Comparison of historical assemblies against the benchmark to include ###
### things like Ash1, HPRCv1, new Verkko, new Hifiasm, and whatever older ###
### HG002 assemblies we can find in GenBank to show progress towards ###
### personalized genomes.] ###
################

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
assemblycolors <- c("#44AA99", "#332288", "#882255", "#888888", "#661100","#999933","#88CCEE", "#CC6677", "#DDCC77")
qvmethodcolors <- c("#DDCC77", "#661100", "#6699CC")
#qvmethodcolors <- c("#661100", "#6699CC")
#assemblynames <- c("ash1v2", "HPRC_year1_mat", "hg002_curated_mat", "hg002_lc24medaka_hap1")
#assemblylabels <- c("Ash1v2.0", "Year 1 HPRC", "HPRC Curated", "LC24 Medaka")
#assemblynames <- c("ash1v2", "HPRC_year1_mat", "hg002_lc24medaka_hap1", "hifi_q28_trio_hic_hap1")
assemblynames <- c("ash1v2", "hprc_year1", "hprc_year2_polished", "hifi_q28_trio_hic")
assemblylabels <- c("Ash1v2 2020", "Yr1 HPRC 2023", "Yr2_polished HPRC 2024", "Revio/Q28 2024")
fiveassemblynames <- c("ash1v2", "hifiasm_2021", "hprc_year1", "hprc_year2_polished", "hifi_q28_trio_hic")
fiveassemblylabels <- c("Ash1v2 2020", "Hifiasm 2021", "Yr1 HPRC 2023", "Yr2 HPRC 2024", "Revio/Q28 2024")
fivewithv2trionames <- c("ash1v2", "hifiasm_2021", "hprc_year1", "hprc_year2_polished", "v2_trio")
fivewithv2triolabels <- c("Ash1v2 2020", "Hifiasm 2021", "Yr1 HPRC 2023", "Yr2 HPRC 2024", "Verkko2 Trio 2025")
ksassemblynames <- c("ash1v2", "hg002v0.1", "hifiasm_2021", "hprc_year1", "hprc_curated", "hprc_year2_polished", "lc24_medaka_6b4_test", "v2_trio", "hprc_year2_v2_thic")
ksassemblylabels <- c("Ash1v2 2020", "HG002v0.1", "Hifiasm 2021", "Yr1 HPRC 2023", "Jarvis HPRC Curated", "HPRC Release2 HiFiasm", "LC24 ONT Medaka", "Verkko2 Trio 2025", "HPRC Release2 Verkko")
adamsassemblynames <- c("ash1v2", "hifiasm_2021", "hprc_year1", "verkko_2023", "lc24_medaka_6b4_test")
adamsassemblylabels <- c("Ash1 2020", "Hifiasm 2021", "HPRCv1 2022", "Verkko 2023", "ONT LC24 2024")
skassemblynames <- c("ash1v2", "verkko_2023", "hifiasm_2021", "hprc_year1", "hprc_year2_polished", "lc24_nopolish", "lc24_medaka_6b4_test", "v2_trio", "hifiasm_ontonly_2025")
skassemblylabels <- c("Ash1v2 2020", "Verkko 2023", "Hifiasm 2021", "Yr1 HPRC 2023", "HPRC Release2 HiFiasm", "LC24 Unpolished", "LC24 ONT Medaka", "Verkko2 Trio 2025", "Hifiasm ONTonly 2025")
test1assemblynames <- c("lc24_nopolish")
test1assemblylabels <- c("LC24 Unpolished")

# NGA plots
assemblysizefiles <- sapply(assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
fiveassemblysizefiles <- sapply(fiveassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
fivewithv2triosizefiles <- sapply(fivewithv2trionames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
ksassemblysizefiles <- sapply(ksassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
adamsassemblysizefiles <- sapply(adamsassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
skassemblysizefiles <- sapply(skassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})

# Make NGAx plots:

pdf("Figure3Output/OriginalFourNGAxPlot.pdf")
assembly_ngax_plot(assemblysizefiles, assemblylabels=assemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
dev.off()

pdf("Figure3Output/FiveAssemblyNGAxPlot.pdf")
assembly_ngax_plot(fiveassemblysizefiles, assemblylabels=fiveassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
dev.off()

pdf("Figure3Output/FiveWithVerkko2TrioNGAxPlot.pdf")
assembly_ngax_plot(fivewithv2triosizefiles, assemblylabels=fiveassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
dev.off()

pdf("Figure3Output/KitchenSinkNGAxPlot.pdf")
assembly_ngax_plot(ksassemblysizefiles, assemblylabels=ksassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
dev.off()

pdf("Figure3Output/AdamsNGAxPlot.pdf")
assembly_ngax_plot(adamsassemblysizefiles, assemblylabels=adamsassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
dev.off()

# Mononucleotide accuracy:
mnstatsfiles <- sapply(assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})
ksmnstatsfiles <- sapply(ksassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})
adamsmnstatsfiles <- sapply(adamsassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})
skmnstatsfiles <- sapply(skassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})

# Make mononucleotide accuracy/QV plots

pdf("Figure3Output/AssemblyMononucAccuracy.pdf")  
assembly_mononucacc_plot(mnstatsfiles, assemblylabels)
dev.off()

pdf("Figure3Output/KitchenSinkMononucAccuracy.pdf")  
assembly_mononucacc_plot(ksmnstatsfiles, ksassemblylabels)
dev.off()

pdf("Figure3Output/KitchenSinkAssemblyMononucErrorQVs.pdf")  
assembly_mononucqv_plot(ksmnstatsfiles, ksassemblylabels)
dev.off()

pdf("Figure3Output/AdamsMononucErrorQVs.pdf")  
assembly_mononucqv_plot(adamsmnstatsfiles, adamsassemblylabels)
dev.off()

# Indel errors:
indelstatsfiles <- sapply(assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})
ksindelstatsfiles <- sapply(ksassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})
adamsindelstatsfiles <- sapply(adamsassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})
skindelstatsfiles <- sapply(skassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})

# Make indels plots:

pdf("Figure3Output/IndelRates.pdf", width=9, height=8)  
assembly_indels_plot(indelstatsfiles, assemblylabels, titleval="Indel discrepancies in assemblies")
dev.off()

pdf("Figure3Output/KitchenSinkIndelRates.pdf", width=9, height=8)  
assembly_indels_plot(ksindelstatsfiles, ksassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=18.0)
dev.off()

pdf("Figure3Output/AdamsRates.pdf", width=9, height=8)  
assembly_indels_plot(adamsindelstatsfiles, adamsassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=18.0)
dev.off()


### Phasing errors plot

phaseswitchfile <- "Figure3Output/switchrates.txt"

pdf("Figure3Output/PhaseSwitchRates.pdf")  
assembly_switchrate_plot(assemblynames, assemblylabels, phaseswitchfile)
dev.off()

pdf("Figure3Output/KitchenSinkPhaseSwitchRates.pdf")  
assembly_switchrate_plot(ksassemblynames, ksassemblylabels, phaseswitchfile)
dev.off()

pdf("Figure3Output/AdamsPhaseSwitchRates.pdf")  
assembly_switchrate_plot(adamsassemblynames, adamsassemblylabels, phaseswitchfile)
dev.off()

### Substitution rates

substitutionstatsfiles <- sapply(assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})
kssubstitutionstatsfiles <- sapply(ksassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})
adamssubstitutionstatsfiles <- sapply(adamsassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})
sksubstitutionstatsfiles <- sapply(skassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})

# Make substitutions plot:

pdf("Figure3Output/SubstitutionRates.pdf", width=9, height=8)  
assembly_substitutions_plot(substitutionstatsfiles, assemblylabels, titleval="Substitution discrepancies in assemblies")
dev.off()

pdf("Figure3Output/KitchenSinkSubstitutionRates.pdf", width=9, height=8)  
assembly_substitutions_plot(kssubstitutionstatsfiles, ksassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
dev.off()

pdf("Figure3Output/AdamsSubstitutionRates.pdf", width=9, height=8)  
assembly_substitutions_plot(adamssubstitutionstatsfiles, adamsassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
dev.off()

### Quality value (QV) scores for different assemblies by different methods
#qvfile <- "Figure3Output/AssemblyQVScores.txt"
yakqvfile <- "Figure3Output/yak_sprq_element_std_hybrid_qvs.txt"
merqqvfile <- "Figure3Output/merqury_qvs.txt"
assemblyqvfile <- "Figure3Output/gqc_assembly_qvs.txt"

pdf("Figure3Output/QVMeasures.pdf", width=13, height=8.5)  
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, assemblynames, assemblylabels)
dev.off()

pdf("Figure3Output/KitchenSinkQVMeasures.pdf", width=13, height=8.5)  
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, ksassemblynames, ksassemblylabels)
dev.off()

pdf("Figure3Output/AdamsQVMeasures.pdf", width=13, height=8.5)  
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, adamsassemblynames, adamsassemblylabels)
dev.off()

### Plot them all together:
pdf("Figure3Multiplot.pdf", width=13, height=9)
par(mfrow=c(2,3))
assembly_ngax_plot(assemblysizefiles, assemblylabels=assemblylabels, plottitle="NGAx for different assemblies")
assembly_qv_plot(qvfile, assemblylabels)
assembly_mononucqv_plot(mnstatsfiles, assemblylabels, pointcex=1.5)
assembly_substitutions_plot(substitutionstatsfiles, assemblylabels, titleval="Substitution discrepancies in assemblies")
assembly_indels_plot(indelstatsfiles, assemblylabels, titleval="Indel discrepancies in assemblies")
assembly_switchrate_plot(phaseswitchfile, plottitle="Phase switch errors per mb")
dev.off()

### Plot them all together:
pdf("Figure3KitchenSinkMultiplot.pdf", width=13, height=9)
par(mfrow=c(2,3))
assembly_ngax_plot(ksassemblysizefiles, assemblylabels=ksassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, ksassemblynames, ksassemblylabels)
assembly_mononucqv_plot(ksmnstatsfiles, ksassemblylabels)
assembly_substitutions_plot(kssubstitutionstatsfiles, ksassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
assembly_indels_plot(ksindelstatsfiles, ksassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=18.0)
assembly_switchrate_plot(ksassemblynames, ksassemblylabels, phaseswitchfile)
dev.off()

### Plot them all together:
pdf("Figure3AdamsMultiplot.pdf", width=13, height=9)
par(mfrow=c(2,3))
assembly_ngax_plot(adamsassemblysizefiles, assemblylabels=adamsassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, adamsassemblynames, adamsassemblylabels)
assembly_mononucqv_plot(adamsmnstatsfiles, adamsassemblylabels, pointcex=1.0)
assembly_substitutions_plot(adamssubstitutionstatsfiles, adamsassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
assembly_indels_plot(adamsindelstatsfiles, adamsassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=20.0)
assembly_switchrate_plot(adamsassemblynames, adamsassemblylabels, phaseswitchfile)
dev.off()

pdf("Figure3AdamsMultiplotLinesBars.pdf", width=13, height=9)
par(mfrow=c(2,3))
assembly_ngax_plot(adamsassemblysizefiles, assemblylabels=adamsassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, adamsassemblynames, adamsassemblylabels)
assembly_mononucqv_plot(adamsmnstatsfiles, adamsassemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
#assembly_substitutions_plot(adamssubstitutionstatsfiles, adamsassemblylabels, legendypos=10, titleval="Substitution discrepancies in assemblies", ybreak=c(40, 120), ymax=160, spanbreak=c(9,12))
assembly_substitutions_plot(adamssubstitutionstatsfiles, adamsassemblylabels, legendypos=10, titleval="Substitution discrepancies in assemblies" )
assembly_indels_plot(adamsindelstatsfiles, adamsassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=20.0)
assembly_switchrate_plot(adamsassemblynames, adamsassemblylabels, phaseswitchfile)
dev.off()

### Plot them all together:
pdf("Figure3SergesKitchenSinkMultiplot.pdf", width=13, height=9)
par(mfrow=c(2,3))
assembly_ngax_plot(skassemblysizefiles, assemblylabels=skassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, skassemblynames, skassemblylabels)
assembly_mononucqv_plot(skmnstatsfiles, skassemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
assembly_substitutions_plot(sksubstitutionstatsfiles, skassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
assembly_indels_plot(skindelstatsfiles, skassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=62.0)
assembly_switchrate_plot(skassemblynames, skassemblylabels, phaseswitchfile)
dev.off()

# see how it looks to plot just one assembly's data:
pdf("Test1AssemblyMultiplot.pdf", width=9, height=9)
par(mfrow=c(2,2))
test1assemblysizefiles <- sapply(test1assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
test1substitutionstatsfiles <- sapply(test1assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})
test1indelstatsfiles <- sapply(test1assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})
test1mnstatsfiles <- sapply(test1assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})
assembly_ngax_plot(test1assemblysizefiles, assemblylabels=test1assemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
assembly_mononucqv_plot(test1mnstatsfiles, test1assemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
assembly_substitutions_plot(test1substitutionstatsfiles, test1assemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
assembly_indels_plot(test1indelstatsfiles, test1assemblylabels, titleval="Indel discrepancies in assemblies", legendypos=62.0)
dev.off()


