##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2021
# Title:       SEE01 paper 1 repeatability plots
# Description: Also include scatterplot of LOD scores vs correlation across all positions for all chemicals, fit a lm, then try doing each chemical separately to see if there are some chemicals where there is a correlation b/w LOD score and per-site correlation. May want to try to weight average correlations across replicates so that haplotypes that moved more are given more weight. Maybe can compute mean correlation across replicates for the top two or three moving haplotypes.

##############################################################################

# ============================================================================
# Load packages and sourced files
# Sourced files are kept in the default working directory	

library(tictoc)
library(data.table)
library(tidyverse)
library(scales)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggbeeswarm)
library(DescTools)
library(abind)
library(job)
library(ggdendro)
library(cowplot)
library(gridGraphics)

source('formatting/Haplotype_file_splitter.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')
source('calculating/Calculating_test_statistics.R')
source('formatting/Haplotype_file_reformatter.R')
source('seq/Position_offsetter.R')



# ============================================================================
# Set global options

defDir <- getwd()
#projectDir <- "/Users/robertlinder/Dropbox/Long_lab/DXQTL03/Primary_experiments/"
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

# ============================================================================
# Custom functions


# ============================================================================
# Load data needed for downstream analyses

founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)
treatKeyDT <- fread("treatment_key.txt", header = T)
offsets <- read.table("newoffsets.txt", header = T)
offsets$lines <- ceiling(offsets[,2]/50)
offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
offsets$chr <- 1:17
g_l <- c(0, cumsum(offsets$len))
ch.bounds <- c(0, g_l[1:17] + offsets[,2])
#correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
max.pos <- g_l[length(g_l)]
mid.ch <- diff(ch.bounds)/2
midpt.ch <- ch.bounds[2:18] - mid.ch
gray <- col2rgb('grey50')
chromBarCol <- rgb(gray[1],gray[2],gray[3], maxColorValue = 255, alpha = 100)
colorCodes = c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
hx = sapply(strsplit(colorCodes, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
hx2 = hx[-c(5,21,24)]   
mycols = c(hx2,"#000000")

# ============================================================================
# Plot a dendrogram of ALL replicates for ALL outbred sexual samples using Spearman correlation; can see from this DMSO R08 doesn't cluster with the rest of the DMSO samples 

directory <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/")
setwd(directory)
allCorDF <- read.table("complete_chems_all_reps_spearman_cors.txt", header = T, sep = "\t")
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

job::job(corJob1 = {
    chemReps <- unique(names(allCorDF))
    corAvgs <- lapply(chemReps, function(x) {
        print(x)
        flush.console()
        tic()
        corrGrp <- allCorDF[,grep(x, names(allCorDF))]
        corrGrp <- as.data.frame(corrGrp)
        names(corrGrp) <- x
        rownames(corrGrp) <- rownames(allCorDF)
        corrGrp$id <- gsub(".*18way_", "", rownames(allCorDF))
        corrGrp$id <- gsub("_12.*-", "-", corrGrp$id)
        #corrGrp$gp <- as.numeric(gsub(".*_12-|-R.*", "", rownames(allCorDF)))
        corrGrpDT <- as.data.table(corrGrp)
        corrAvgs <- corrGrpDT[, mean(get(x), na.rm = T), by = id]
        setnames(corrAvgs, c("id", x))
        toc()
        corrAvgs
    } )
}, import = c(allCorDF), packages = c("data.table", "tictoc") )

chems <- fread("/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/Data_Storage/Sequencing_Data_Processed/chemAbbrev.txt", sep = "\t", header = T)

sSamples <- mrgeData[Type == "Outbred_sexual"]
samplesChems <- merge(sSamples, chems, by.x = "chem", by.y = "abbr", all.x = T)
samplesChems[, replicate := substr(id, 14, 16)]
samplesChems[, newId := paste0(chemical, "-", replicate)]
samplesChems[, useId := paste0(chem, "-", replicate)]

filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/")
CC2 <- do.call(cbind, corJob1$corAvgs)
CC3 <- as.data.frame(CC2)
rownames(CC3) <- CC2$id
CC4 <- CC3[, -grep("^id", names(CC3))]
names(CC4) <- gsub("X18way_", "", names(CC4))
names(CC4) <- gsub("_12.", "-", names(CC4))
CC5 <- CC4[, names(CC4) %in% samplesChems$newId]
CC6 <- CC5[rownames(CC5) %in% samplesChems$newId,]
nameDT <- data.table(names = names(CC6), ids = samplesChems$newId) 
checkOrder <- all(nameDT$names == nameDT$ids)
names(CC6) <- samplesChems$useId

rowNameDT <- data.table(names = rownames(CC6), ids = samplesChems$newId) 
checkOrder <- all(rowNameDT$names == rowNameDT$ids)
rownames(CC6) <- samplesChems$useId

out <- hclust(as.dist(1-CC6^2))

fileName <- paste0(filePath, "complete_sexual_reps_all_chems_Spearman_cors_revised.pdf")
mainTitle <- paste0("Spearman_correlations")
pdf(fileName, width = 16, height = 10)
par(mar = c(3,3,1.5,0.5), mgp = c(2,0.7,0), oma = c(0,0,0,0))
par(cex = 1)
plot(as.dendrogram(out), xlab = "", ylab = "", main = "", sub = "", axes = FALSE, ylim = c(0, 1))
#rect.hclust(out, h = 0.85, border = "red")
par(cex = 1)
title(main = mainTitle, ylab = "height")
axis(2)
dev.off()


# ============================================================================
# Plot a dendrogram of ALL replicates for ALL sample using Spearman correlation; exclude chlorpromazine R5 (low heterozygosity)

directory <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/")
setwd(directory)
allCorDF <- read.table("complete_chems_all_reps_spearman_cors.txt", header = T, sep = "\t")
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

job::job(corJob1 = {
chemReps <- unique(names(allCorDF))
corAvgs <- lapply(chemReps, function(x) {
    print(x)
    flush.console()
    tic()
    corrGrp <- allCorDF[,grep(x, names(allCorDF))]
    corrGrp <- as.data.frame(corrGrp)
    names(corrGrp) <- x
    rownames(corrGrp) <- rownames(allCorDF)
    corrGrp$id <- gsub(".*18way_", "", rownames(allCorDF))
    corrGrp$id <- gsub("_12.*-", "-", corrGrp$id)
    #corrGrp$gp <- as.numeric(gsub(".*_12-|-R.*", "", rownames(allCorDF)))
    corrGrpDT <- as.data.table(corrGrp)
    corrAvgs <- corrGrpDT[, mean(get(x), na.rm = T), by = id]
    setnames(corrAvgs, c("id", x))
    toc()
    corrAvgs
} )
}, import = c(allCorDF), packages = c("data.table", "tictoc") )


filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/")
CC2 <- do.call(cbind, corJob1$corAvgs)
CC3 <- as.data.frame(CC2)
rownames(CC3) <- CC2$id
CC4 <- CC3[, -grep("^id", names(CC3))]
names(CC4) <- gsub("X18way_", "", names(CC4))
names(CC4) <- gsub("_12.", "-", names(CC4))
CC5 <- CC4[, -grep("sulfite|tunicamycin|caffeine|cisplatin|ethanol|dmso|sodium_chloride-R13|sodium_chloride-R16|chlorpromazine-R05", names(CC4))]
CC6 <- CC5[-grep("sulfite|tunicamycin|caffeine|cisplatin|ethanol|dmso|sodium_chloride-R13|sodium_chloride-R16|chlorpromazine-R05", rownames(CC5)),]

CC7 <- CC6[-grep("fluconazole|nic|cadmium_chloride-R0[1-3]|YPD-R10|glacial_acetic_acid-R08|chlorpromazine-R12|chlorpromazine-R13", rownames(CC6)),]
CC8 <- CC7[,-grep("fluconazole|nic|cadmium_chloride-R0[1-3]|YPD-R10|glacial_acetic_acid-R08|chlorpromazine-R12|chlorpromazine-R13", names(CC7))]

write.table(CC6, "sexual_reps_6_plus_spearman_avgs.txt", quote = F, sep = "\t", col.names = T)

CC6 <- read.table("sexual_reps_6_plus_spearman_avgs.txt", header = T, sep = "\t", as.is = T)
out <- hclust(as.dist(1-CC6^2))
#outSelect55 <- hclust(as.dist(1-CC8^2))

fileName <- paste0(filePath, "all_reps_dendrogram_revised.pdf")
#mainTitle <- paste0("Spearman_correlations")
pdf(fileName, width = 8, height = 5)
par(mar = c(8.7,3,1.5,0.5), mgp = c(2,0.7,0), oma = c(0,0,0,0))
par(cex = 0.75)
plot(as.dendrogram(out), xlab = "", ylab = "", main = "", sub = "", axes = FALSE, ylim = c(0, 1))
#rect.hclust(out, h = 0.85, border = "red")
par(cex = 1)
#title(main = mainTitle, ylab = "height")
axis(2)
allrepplot <- recordPlot()
#plot.new()
dev.off()

fileName2 <- paste0(filePath, "selected_reps_dendrogram_revised.pdf")

pdf(fileName2, width = 8, height = 5)
par(mar = c(8.7,3,1.5,0.5), mgp = c(2,0.7,0), oma = c(0,0,0,0))
par(cex = 0.75, xpd = FALSE)
plot(as.dendrogram(outSelect55), xlab = "", ylab = "", main = "", sub = "", axes = FALSE, ylim = c(0, 1))
rect.hclust(outSelect55, h = 0.85, border = "red")
par(cex = 1)
#title(main = mainTitle, ylab = "height")
axis(2)
somerepplot <- recordPlot()
#plot.new()
dev.off()

#savePlot <- plot_grid(allrepplot, somerepplot,labels = c("A", "B"), ncol = 1, align = "v", rel_heights = c(2,1))
#ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS5_dendrograms.png"), savePlot, width = 8, height = 10, units = "in", dpi = 350)


# ============================================================================
# Plot a dendrogram of selected replicates using Spearman correlation

directory <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/")
setwd(directory)
allCorDF <- read.table("complete_chems_all_reps_spearman_cors.txt", header = T, sep = "\t")
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

allCorDF2 <- allCorDF[-grep("fluconazole|nic", rownames(allCorDF)),]
allCorDF2 <- allCorDF2[-grep("cadmium_chloride_12.*-R0[1-3]|YPD_12.*-R10|glacial_acetic_acid_12.*-R08|chlorpromazine_12.*-R12|chlorpromazine_12.*-R13", rownames(allCorDF2)),]
allCorDF2 <- allCorDF2[, -grep("fluconazole|nic|cadmium_chloride_12.R0[1-3]|YPD_12.R10|glacial_acetic_acid_12.R08|chlorpromazine_12.R12|chlorpromazine_12.R13", names(allCorDF2))]

job::job(corJob = {
    chemReps <- unique(names(allCorDF2))
    corAvgs <- lapply(chemReps, function(x) {
        tic()
        corrGrp <- allCorDF2[,grep(x, names(allCorDF2))]
        corrGrp <- as.data.frame(corrGrp)
        names(corrGrp) <- x
        rownames(corrGrp) <- rownames(allCorDF2)
        corrGrp$id <- gsub(".*18way_", "", rownames(allCorDF2))
        corrGrp$id <- gsub("_12.*-", "-", corrGrp$id)
        #corrGrp$gp <- as.numeric(gsub(".*_12-|-R.*", "", rownames(allCorDF)))
        corrGrpDT <- as.data.table(corrGrp)
        corrAvgs <- corrGrpDT[, mean(get(x), na.rm = T), by = id]
        setnames(corrAvgs, c("id", x))
        toc()
        corrAvgs
    } )
}, import = c(allCorDF2), packages = c("data.table", "tictoc") )


# filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/")
# CC2 <- do.call(cbind, corJob$corAvgs)
# CC3 <- as.data.frame(CC2)
# rownames(CC3) <- CC2$id
# CC4 <- CC3[, -grep("^id", names(CC3))]
# names(CC4) <- gsub("X18way_", "", names(CC4))
# names(CC4) <- gsub("_12.", "-", names(CC4))
# out <- hclust(as.dist(1-CC4^2))

fileName <- paste0(filePath, "selected_reps_dendrogram.pdf")

pdf(fileName, width = 8, height = 5)
par(mar = c(8.7,3,1.5,0.5), mgp = c(2,0.7,0), oma = c(0,0,0,0))
par(cex = 0.75)
plot(as.dendrogram(out), xlab = "", ylab = "", main = "", sub = "", axes = FALSE, ylim = c(0, 1))
rect.hclust(out, h = 0.85, border = "red")
par(cex = 1)
#title(main = mainTitle, ylab = "height")
axis(2)
#somerepplot <- recordPlot()
dev.off()

#savePlot <- plot_grid(allrepplot, somerepplot,labels = c("A", "B"), ncol = 1, hjust = -1)

# ============================================================================
# Trouble-shooting