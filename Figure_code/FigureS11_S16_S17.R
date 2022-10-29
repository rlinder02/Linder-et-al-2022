##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2019
# Title:       Plotting raw haplotype frequencies 
# Description: Plotting haplotype frequencies of the base population, along with the treatment replicates and ypd-only replicates; also making a table of raw haplotype frequencies and heterozygosity 

##############################################################################

# ============================================================================
# Load packages and sourced files
# Sourced files are kept in the default working directory
library(tictoc)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(readxl)
library(gprofiler2)
library(gridExtra)
library(egg)
library(scales)
library(dplyr)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggcorrplot)
library(ggbeeswarm)
library(tidyverse)
library(DescTools)
library(lemon)
library(abind)
library(cowplot)
library(tictoc)
library(hexbin)
library(ggpmisc)

source('formatting/Haplotype_file_splitter.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')
source('calculating/Calculating_test_statistics.R')

# ============================================================================
# Set global options

defDir <- getwd()
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))
options("scipen"=999, "digits" = 4)
# ============================================================================
# Custom functions

add.alpha <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
              rgb(x[1], x[2], x[3], alpha=alpha)) 
}    
rep.list <- function(object, repObject) {
    rep(list(object), length(repObject))
}

# ============================================================================
# Load data

avgHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Averaged_and_sd_tx_tables/", analysisType = "hap_freqs_avg_sd_difs", samplePattern = ".*")
avgHapDifsDTsCpy <- copy(avgHapDifsDTs)

# ============================================================================
# Find the two most increased haplotypes per position.

idxLooper <- lapply(avgHapDifsDTs, function(chem)  {
    freqDifs <- chem[, tstrsplit(avgDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
    baseFreqs <- chem[, tstrsplit(baseFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
    avgFreqs <- as.data.frame(chem[, tstrsplit(avgFreq, split = ";", type.convert = TRUE, fixed = TRUE)])
    collHaps <- as.data.frame(chem[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
    baseFreqsDF <- as.data.frame(baseFreqs)
    maxVals <- apply(freqDifs, 1, max, na.rm = T) ## only looking at increasing values
    maxIdx <- apply(freqDifs, 1, which.max)
    maxHap <- collHaps[cbind(seq_along(maxIdx), maxIdx)]
    baseHap <- baseFreqsDF[cbind(seq_along(maxIdx), maxIdx)]
    evFreq <- avgFreqs[cbind(seq_along(maxIdx), maxIdx)]
    chem[, c("maxHap", "maxChange", "baseHap", "evFreq") := .(maxHap, maxVals, baseHap, evFreq)]
    nxtMax <- apply(freqDifs, 1, function(x) rev(sort(x))[2])
    nxtMaxIdx <- apply(freqDifs, 1, function(x) which(x == rev(sort(x))[2])[1])
    maxHap2 <- collHaps[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    baseHap2 <- baseFreqsDF[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    evFreq2 <- avgFreqs[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    chem[, c("maxHap2", "maxChange2", "baseHap2", "evFreq2") := .(maxHap2, nxtMax, baseHap2, evFreq2)]
    idx <- chem[, c("chr", "pos", "gp", "Chemical", "maxHap", "maxChange", "maxHap2", "maxChange2", "baseHap", "baseHap2", "evFreq", "evFreq2")]
    idx[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
    idxMelt1 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("maxChange", "maxChange2"))
    idxMelt2 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("maxHap", "maxHap2"))
    idxMelt3 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("baseHap", "baseHap2"), variable.name = "maxBase", value.name =  "baseHapFreq")
    idxMelt4 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("evFreq", "evFreq2"), variable.name = "maxEv", value.name =  "evFreq")
    setnames(idxMelt2, old = c("variable", "value"), new = c("hapType", "haplotype"))
    idxMrge <- idxMelt1[, c("hapType", "haplotype", "maxBase", "baseHapFreq", "maxEv", "evFreq") := .(idxMelt2$hapType, idxMelt2$haplotype, idxMelt3$maxBase, idxMelt3$baseHapFreq, idxMelt4$maxEv, idxMelt4$evFreq)]
    setnames(idxMrge, old = c("variable", "value"), new = c("changeType", "frequencyChangeNorm"))
    idxMrge
})

chemMrge <- do.call(rbind, idxLooper)
setnames(chemMrge, old = c("frequencyChangeNorm", "haplotype", "baseHapFreq"), new = c("frequencyChange", "haplotype2", "baseHapFreq2"))
baseFreqCnts <- chemMrge[,.(count = .N), by = baseHapFreq2][, percent := prop.table(count)*100][order(-count),][]

chemDT <- dcast(chemMrge, gp + chr + pos + Chemical ~ changeType, value.var = "frequencyChange")
dummy <- chemDT[, .(maxChange = max(maxChange)), by = Chemical][, maxChange2 := maxChange]

avgMaxChngeGW <- chemDT[, mean(maxChange)]
avgMaxChnge2GW <- chemDT[, mean(maxChange2)]

maxToNxtMax <- avgMaxChngeGW/avgMaxChnge2GW

# ============================================================================
# Find the two most increased haplotypes at positions at which the second most changed haplotype started at a higher frequency than the most changed haplotype

idxLooper2 <- lapply(avgHapDifsDTsCpy, function(chem)  {
    freqDifs <- chem[, tstrsplit(avgDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
    baseFreqs <- chem[, tstrsplit(baseFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
    collHaps <- as.data.frame(chem[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
    baseFreqsDF <- as.data.frame(baseFreqs)
    basePQ <- baseFreqs *(1-baseFreqs)
    freqDifsNorm <- freqDifs/basePQ
    freqDifsNormDF <- as.data.frame(freqDifsNorm)
    maxVals <- apply(freqDifs, 1, max, na.rm = T) ## only looking at increasing values
    maxIdx <- apply(freqDifs, 1, which.max)
    maxHap <- collHaps[cbind(seq_along(maxIdx), maxIdx)]
    baseHap <- baseFreqsDF[cbind(seq_along(maxIdx), maxIdx)]
    maxS <- freqDifsNormDF[cbind(seq_along(maxIdx), maxIdx)]
    chem[, c("maxHap", "maxChange", "maxS", "baseHap") := .(maxHap, maxVals, maxS, baseHap)]
    nxtMax <- apply(freqDifs, 1, function(x) rev(sort(x))[2])
    nxtMaxIdx <- apply(freqDifs, 1, function(x) which(x == rev(sort(x))[2])[1])
    maxHap2 <- collHaps[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    baseHap2 <- baseFreqsDF[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    maxS2 <- freqDifsNormDF[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    chem[, c("maxHap2", "maxChange2", "maxS2", "baseHap2") := .(maxHap2, nxtMax, maxS2, baseHap2)]
    idx <- chem[, c("chr", "pos", "gp", "Chemical", "maxHap", "maxChange", "maxHap2", "maxChange2", "baseHap", "baseHap2", "maxS", "maxS2")]
    idx[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
    idx <- idx[baseHap2 > baseHap]
    idxMelt1 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("maxChange", "maxChange2"))
    idxMelt2 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("maxHap", "maxHap2"))
    idxMelt3 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("baseHap", "baseHap2"), variable.name = "maxBase", value.name =  "baseHapFreq")
    idxMelt4 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("maxS", "maxS2"), variable.name = "hapS", value.name =  "S")
    setnames(idxMelt2, old = c("variable", "value"), new = c("hapType", "haplotype"))
    idxMrge <- idxMelt1[, c("hapType", "haplotype", "maxBase", "baseHapFreq", "hapS", "S") := .(idxMelt2$hapType, idxMelt2$haplotype, idxMelt3$maxBase, idxMelt3$baseHapFreq, idxMelt4$hapS, idxMelt4$S)]
    setnames(idxMrge, old = c("variable", "value"), new = c("changeType", "frequencyChangeNorm"))
    idxMrge
})

chemMrge2 <- do.call(rbind, idxLooper2)
setnames(chemMrge2, old = c("frequencyChangeNorm", "haplotype", "baseHapFreq"), new = c("frequencyChange3", "haplotype3", "baseHapFreq3"))
baseFreqCnts2 <- chemMrge2[,.(count = .N), by = baseHapFreq3][, percent := prop.table(count)*100][order(-count),][]

chemDT2 <- dcast(chemMrge2, gp + chr + pos + Chemical ~ changeType, value.var = "frequencyChange3")
dummy2 <- chemDT2[, .(maxChange = max(maxChange)), by = Chemical][, maxChange2 := maxChange]

avgMaxChngeGW2 <- chemDT2[, mean(maxChange)]
avgMaxChnge2GW2 <- chemDT2[, mean(maxChange2)]

maxToNxtMax2 <- avgMaxChngeGW2/avgMaxChnge2GW2

# ============================================================================
# Plot most changed vs next most changed haplotype frequencies per locus (A) and at loci where the next most changed haplotype starts at a higher frequency than the most changed haplotype (B)

panelPlotA <- ggplot(chemDT, aes(maxChange, maxChange2)) + xlab("Average change of most increased haplotype") + ylab("Average change of 2nd most increased haplotype") + facet_wrap(~Chemical, scales = "free", nrow = 2) + geom_blank(data = dummy) + geom_abline(intercept = 0, slope = 1, colour = "black", lty = 3) + geom_smooth(span = 0.3, colour = "red") + geom_hex(bins = 50, alpha = 0.75) + scale_fill_continuous(type = "viridis") + theme_bw(base_size = 10) + theme(panel.grid = element_blank(), aspect.ratio = 1)  + scale_x_continuous(limits = c(0, NA), expand = c(0,0)) + scale_y_continuous(limits = c(0, NA), expand = c(0,0))

panelPlotB <- ggplot(chemDT2, aes(maxChange, maxChange2)) + xlab("Average change of most increased haplotype") + ylab("Average change of 2nd most increased haplotype") + facet_wrap(~Chemical, scales = "free", nrow = 2) + geom_blank(data = dummy2) + geom_abline(intercept = 0, slope = 1, colour = "black", lty = 3) + geom_smooth(span = 0.3, colour = "red") + geom_hex(bins = 50, alpha = 0.75) + scale_fill_continuous(type = "viridis") + theme_bw(base_size = 10) + theme(panel.grid = element_blank(), aspect.ratio = 1)  + scale_x_continuous(limits = c(0, NA), expand = c(0,0)) + scale_y_continuous(limits = c(0, NA), expand = c(0,0))

savePlot <- plot_grid(panelPlotA, panelPlotB, align = 'vh', labels = c("A", "B"), ncol = 1, hjust = -1)

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS11_1st_2nd_increased_haps_gwide.pdf"), savePlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS11_1st_2nd_increased_haps_gwide.png"), savePlot, width = 8, height = 10, units = "in", dpi = 350)

system2('pdfcrop', c('filename1', 'filename2'))
# ============================================================================
# Plot the frequency in the base vs the average over replicates of the most changed haplotype (A) and the next-most changed haplotype (B) genomewide

panelPlotA <- ggplot(chemMrge[changeType == "maxChange"], aes(baseHapFreq2, evFreq)) + xlab("Base frequency") + ylab("Average evolved frequency") + facet_wrap(~Chemical, scales = "free_x", nrow = 2) + geom_smooth(method = "lm", se=FALSE, colour = "red", size = 0.5) + stat_regline_equation(aes(label = ..rr.label..), label.x.npc = "right", label.y.npc = "bottom", hjust = 1, vjust = 0.2, size = 3) + geom_hex(bins = 50, alpha = 0.75) + scale_fill_continuous(type = "viridis") + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + coord_cartesian(ylim = c(0, 1))  

panelPlotB <- ggplot(chemMrge[changeType == "maxChange2"], aes(baseHapFreq2, evFreq)) + xlab("Base frequency") + ylab("Average evolved frequency") + facet_wrap(~Chemical, scales = "free_x", nrow = 2) + geom_smooth(method = "lm", se=FALSE, colour = "red", size = 0.5) + stat_regline_equation(aes(label = ..rr.label..), label.x.npc = "right", label.y.npc = "bottom", hjust = 1, vjust= 0.2, size = 3) + geom_hex(bins = 50, alpha = 0.75) + scale_fill_continuous(type = "viridis") + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + coord_cartesian(ylim = c(0, 1)) 

savePlot <- plot_grid(panelPlotA, panelPlotB, align = 'vh', labels = c("A", "B"), ncol = 1)

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS16_1st_2nd_increased_haps_freqs_top_peaks.pdf"), savePlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS16_1st_2nd_increased_haps_freqs_top_peaks.png"), savePlot, width = 8, height = 10, units = "in", dpi = 350)

# ============================================================================
# Plot the frequency in the base vs the average change over replicates of the most changed haplotype (A) and the next-most changed haplotype (B) genomewide

my.formula <- y ~ x

panelPlotA <- ggplot(chemMrge[changeType == "maxChange"], aes(baseHapFreq2, frequencyChange)) + xlab("Base frequency") + ylab("Average haplotype change") + facet_wrap(~Chemical, scales = "free_x", nrow = 2) + geom_smooth(method = "lm", se=FALSE, colour = "red", size = 0.5) + stat_poly_eq(formula = my.formula, aes(label = ..rr.label..), parse = TRUE, label.x = 1, label.y = 1, size = 3, rr.digits = 2) + geom_hex(bins = 50, alpha = 0.75) + scale_fill_continuous(type = "viridis") + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + coord_cartesian(ylim = c(0, 1), xlim = c(0, 0.7))  

panelPlotB <- ggplot(chemMrge[changeType == "maxChange2"], aes(baseHapFreq2, frequencyChange)) + xlab("Base frequency") + ylab("Average haplotype change") + facet_wrap(~Chemical, scales = "free_x", nrow = 2) + geom_smooth(method = "lm", se=FALSE, colour = "red", size = 0.5) + stat_poly_eq(formula = my.formula, aes(label = ..rr.label..), parse = TRUE, label.x = 1, label.y = 1, size = 3, rr.digits = 2) + geom_hex(bins = 50, alpha = 0.75) + scale_fill_continuous(type = "viridis") + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + coord_cartesian(ylim = c(0, 1), xlim = c(0, 0.7)) 

savePlot <- plot_grid(panelPlotA, panelPlotB, align = 'vh', labels = c("A", "B"), ncol = 1)

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS17_1st_2nd_increased_haps_basevschange_gwide.pdf"), savePlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS17_1st_2nd_increased_haps_basevschange_gwide.png"), savePlot, width = 8, height = 10, units = "in", dpi = 350)

# ============================================================================
# Trouble-shooting