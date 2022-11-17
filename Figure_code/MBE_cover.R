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
# Load data needed for downstream analyses

founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)
indHaps <- read.table("Individual_haplotype_differences_55_selected_populations.txt",header=TRUE)
indHapsDT <- as.data.table(indHaps)
avgHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Averaged_and_sd_tx_tables/", analysisType = "hap_freqs_avg_sd_difs", samplePattern = ".*")[[1]]
avgHapDifsDTsCpy <- copy(avgHapDifsDTs)


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
treatmentAbbrv <- fread("chemAbbrev.txt")
plottingFactors <- read.table("all_chems_list.txt", header = T, sep = "\t")
plottingFactors2 <- plottingFactors$Chemical
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
colPal <- c("dodgerblue2", "blue2", "skyblue3", "darkturquoise",
            "brown", "olivedrab", "tomato",
            "chocolate", "burlywood4",
            "saddlebrown",
            "magenta4", "maroon2",
            "mediumpurple",
            "darkcyan",
            "cornsilk4",
            "mediumvioletred")

names(colPal) <- plottingFactors
treatmentAbbrv <- fread("chemAbbrev.txt")
names(colPal) <- treatmentAbbrv$abbr
colorPal <- colPal[order(names(colPal))]
plottingFactors <- read.table("all_chems_list.txt", header = T, sep = "\t")
plottingFactors2 <- plottingFactors$Chemical
plottingFactors3 <- plottingFactors2[c(2,1,4,3,5,6,7,8,9,12,10,11,13,14,15,16)]
names(colorPal) <- plottingFactors3

plottingColors <- colPal
colorDF <- data.frame(factors = plottingFactors2, colors = plottingColors, stringsAsFactors = FALSE)
myColors <- colorDF$colors
names(colPal) <- plottingFactors2
colorCodes = c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
hx = sapply(strsplit(colorCodes, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
hx2 = hx[-c(5,21,24)]   
mycols = c(hx2,"#000000")
mypal <- RColorBrewer::brewer.pal(4, "Dark2")
myColPal <- add.alpha(mycols, 0.5)
topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))

# ============================================================================
# Find the two most increased haplotypes per position.

chem <- indHapsDT
#freqDifs <- chem[, tstrsplit(freqDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
baseFreqs <- as.data.frame(avgHapDifsDTs[, tstrsplit(baseFreqs, split = ";", type.convert = TRUE, fixed = TRUE)])
baseHaps <- as.data.frame(avgHapDifsDTs[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
maxBaseVals <- apply(baseFreqs, 1, max, na.rm = T)
maxBaseIdx <- apply(baseFreqs, 1, which.max)
maxBaseHap <- baseHaps[cbind(seq_along(maxBaseIdx), maxBaseIdx)]
maxBaseFreqs <- baseFreqs[cbind(seq_along(maxBaseIdx), maxBaseIdx)]
avgHapDifsDTs[, c("maxHap", "maxFreq") := .(maxBaseHap, maxBaseFreqs)]
baseDT <- avgHapDifsDTs[, c("chr", "pos", "Chemical", "Replicate", "maxHap", "maxFreq")]
baseDT[, c("Chemical", "Replicate") := .("Base_population", "Base_population")]

evFreqs <- as.data.frame(chem[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)])
evHaps <- as.data.frame(chem[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
evMaxVals <- apply(evFreqs, 1, max, na.rm = T) ## only looking at increasing values
evMaxIdx <- apply(evFreqs, 1, which.max)
evMaxHap <- evHaps[cbind(seq_along(evMaxIdx), evMaxIdx)]
evMaxFreqs <- evFreqs[cbind(seq_along(evMaxIdx), evMaxIdx)]
chem[, c("maxHap", "maxFreq") := .(evMaxHap, evMaxFreqs)]
evDT <- chem[, c("chr", "pos", "Chemical", "Replicate", "maxHap", "maxFreq")]

mrgeDT <- rbind(baseDT, evDT)

hapColors <- match(mrgeDT$maxHap, topHaps)
hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
Colors <- mycols[hapColors]
mrgeDT[, "colorCodes":= Colors][, "pos(kb)" := pos/1000]
mrgeDT[, Chemical := gsub("18way_", "", Chemical)]
fwrite(mrgeDT, "Most_frequent_haplotype.txt", sep = "\t", col.names = T)

# ============================================================================
# Plot the most frequent haplotype at each position at a specific chromosome with y axes; chromosome 9


mrgeDT <- fread("/Users/robertlinder/Dropbox/Long_lab/SEE01/All_figures/MBE_cover/Most_frequent_haplotype.txt.gz", header = T, sep = '\t')

ksmoother <- function(DT) {
    smooth <- ksmooth(DT$pos, DT$maxFreq, kernel = "normal", bandwidth = 5000)
    DT2 <- copy(DT)
    DT2[, c("pos", "maxFreq") := .(smooth$x, smooth$y)]
    return(DT2)
}

mrgeDTfiltered <- mrgeDT[colorCodes != "#000000"]

chrom <- 9
chrDT <- mrgeDTfiltered[chr == chrom]
baseDT <- chrDT[Chemical == "Base_population"]
cadR1DT <- chrDT[Chemical == "cadmium_chloride" & Replicate == "R04"]
cadRepsDT <- chrDT[Chemical == "cadmium_chloride" & Replicate %in% c("R05", "R10", "R11")]

#cadRepsDT <- chrDT[Chemical == "cadmium_chloride"]

chlorpR1DT <- chrDT[Chemical == "chlorpromazine" & Replicate == "R01"]
chlorpRepsDT <- chrDT[Chemical == "chlorpromazine" & Replicate %in% c("R02", "R04", "R06")]
#chlorpRepsDT <- chrDT[Chemical == "chlorpromazine"]

baseDTsmooth <- ksmoother(baseDT)
cadR1DTsmooth <- ksmoother(cadR1DT)
chlorpR1DTsmooth <- ksmoother(chlorpR1DT)

filePath <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/All_figures/MBE_cover/"
textSize <- 1.5
ymax <- 1
ymin <- 0

pdf(paste0(filePath, "MBE_cover_submission_v1.pdf"), width = 8, height = 10)
#png(paste0(filePath, "MBE_cover_submission_v1.png"), width = 8, height = 10, units = 'in', res = 350)
layout(matrix(1:5))
par(mar = c(0.7,2.5,0,6), mgp = c(1.5,0.7,0), oma = c(2,2.5,2,2))
# plotting base population

baseDTsmooth[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
baseDTsmooth[, points(pos, maxFreq, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
baseDTsmooth[, lines(pos, maxFreq, col = "lightgray", bty = "n", lwd = 0.3)]
baseDTsmooth[, mtext(text="Ancestral population", side = 3, cex = textSize, xpd = NA, line = -1)]
baseDTsmooth[, axis(2, at=c(0,1), labels=c(0,1), las = 1, cex.axis = 0.8)]
# plotting one cadmium chloride replicate population

cadR1DTsmooth[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
cadR1DTsmooth[, points(pos, maxFreq, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
cadR1DTsmooth[, lines(pos, maxFreq, col = "lightgray", bty = "n", lwd = 0.3)]
cadR1DTsmooth[, mtext(text="Cadmium chloride evolved", side = 3, cex = textSize, xpd = NA, line = -1)]
cadR1DTsmooth[, axis(2, at=c(0,1), labels=c(0,1), las = 1, cex.axis = 0.8)]
# plotting three additional cadmium chloride replicate populations

cadRepsDT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
counter <- 0
lapply(unique(cadRepsDT$Replicate), function(p) {
    cadRepsDTsmooth <- ksmoother(cadRepsDT[Replicate == p])
    cadRepsDTsmooth[Replicate == p, points(pos, maxFreq+counter, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
    cadRepsDTsmooth[Replicate == p, lines(pos, maxFreq+counter, col = "lightgray", bty = "n", lwd = 0.3)]
    counter <<- counter + 0.2
} )
# plotting one chlorpromazine replicate population

chlorpR1DTsmooth[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
chlorpR1DTsmooth[, points(pos, maxFreq, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
chlorpR1DTsmooth[, lines(pos, maxFreq, col = "lightgray", bty = "n", lwd = 0.3)]
chlorpR1DTsmooth[, mtext(text="Chlorpromazine evolved", side = 3, cex = textSize, xpd = NA, line = -1)]
chlorpR1DTsmooth[, axis(2, at=c(0,1), labels=c(0,1), las = 1, cex.axis = 0.8)]
# plotting three additional chlorpromazine replicate populations

chlorpRepsDT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
counter <- 0
lapply(unique(chlorpRepsDT$Replicate), function(p) {
    chlorpRepsDTsmooth <- ksmoother(chlorpRepsDT[Replicate == p])
    chlorpRepsDTsmooth[Replicate == p, points(pos, maxFreq+counter, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
    chlorpRepsDTsmooth[Replicate == p, lines(pos, maxFreq+counter, col = "lightgray", bty = "n", lwd = 0.3)]
    counter <<- counter + 0.2
} )
chlorpRepsDTsmooth[, mtext(text=paste0("Chr", as.roman(chrom)), side = 1, line = 1, cex = textSize, xpd = NA)]
dev.off()

# ============================================================================
# Plot the most frequent haplotype at each position at a specific chromosome with y axes; chromosome 16


mrgeDT <- fread("/Users/robertlinder/Dropbox/Long_lab/SEE01/All_figures/MBE_cover/Most_frequent_haplotype.txt.gz", header = T, sep = '\t')

ksmoother <- function(DT) {
    smooth <- ksmooth(DT$pos, DT$maxFreq, kernel = "normal", bandwidth = 5000)
    DT2 <- copy(DT)
    DT2[, c("pos", "maxFreq") := .(smooth$x, smooth$y)]
    return(DT2)
}

mrgeDTfiltered <- mrgeDT[colorCodes != "#000000"]

chrom <- 16
chrDT <- mrgeDTfiltered[chr == chrom]
baseDT <- chrDT[Chemical == "Base_population"]
cadR1DT <- chrDT[Chemical == "cadmium_chloride" & Replicate == "R04"]
cadRepsDT <- chrDT[Chemical == "cadmium_chloride" & Replicate %in% c("R05", "R10", "R11")]

chlorpR1DT <- chrDT[Chemical == "chlorpromazine" & Replicate == "R01"]

chlorpRepsDT <- chrDT[Chemical == "chlorpromazine" & Replicate %in% c("R02", "R04", "R06")]
#chlorpRepsDT <- chrDT[Chemical == "chlorpromazine"]

baseDTsmooth <- ksmoother(baseDT)
cadR1DTsmooth <- ksmoother(cadR1DT)
chlorpR1DTsmooth <- ksmoother(chlorpR1DT)

filePath <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/All_figures/MBE_cover/"
textSize <- 1.5
ymax <- 1
ymin <- 0

#pdf(paste0(filePath, "MBE_cover_submission_v1.pdf"), width = 8, height = 10)
png(paste0(filePath, "MBE_cover_submission_v2.png"), width = 8, height = 10, units = 'in', res = 350)
layout(matrix(1:5))
par(mar = c(0.7,2.5,0,6), mgp = c(1.5,0.7,0), oma = c(2,2.5,2,2))
# plotting base population

baseDTsmooth[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
baseDTsmooth[, points(pos, maxFreq, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
baseDTsmooth[, lines(pos, maxFreq, col = "lightgray", bty = "n", lwd = 0.3)]
baseDTsmooth[, mtext(text="Ancestral population", side = 3, cex = textSize, xpd = NA, line = -1)]
baseDTsmooth[, axis(2, at=c(0,1), labels=c(0,1), las = 1, cex.axis = 0.8)]
# plotting one cadmium chloride replicate population

cadR1DTsmooth[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
cadR1DTsmooth[, points(pos, maxFreq, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
cadR1DTsmooth[, lines(pos, maxFreq, col = "lightgray", bty = "n", lwd = 0.3)]
cadR1DTsmooth[, mtext(text="Cadmium chloride evolved", side = 3, cex = textSize, xpd = NA, line = -1)]
cadR1DTsmooth[, axis(2, at=c(0,1), labels=c(0,1), las = 1, cex.axis = 0.8)]
# plotting two additional cadmium chloride replicate populations
## "R04" "R05" "R06" "R08" "R09" "R10" "R11"

cadRepsDT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
counter <- 0
lapply(unique(cadRepsDT$Replicate), function(p) {
    cadRepsDTsmooth <- ksmoother(cadRepsDT[Replicate == p])
    cadRepsDTsmooth[Replicate == p, points(pos, maxFreq-counter, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
    cadRepsDTsmooth[Replicate == p, lines(pos, maxFreq-counter, col = "lightgray", bty = "n", lwd = 0.3)]
    counter <<- counter + 0.3
} )
# plotting one chlorpromazine replicate population

chlorpR1DTsmooth[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
chlorpR1DTsmooth[, points(pos, maxFreq, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
chlorpR1DTsmooth[, lines(pos, maxFreq, col = "lightgray", bty = "n", lwd = 0.3)]
chlorpR1DTsmooth[, mtext(text="Chlorpromazine evolved", side = 3, cex = textSize, xpd = NA, line = -1)]
chlorpR1DTsmooth[, axis(2, at=c(0,1), labels=c(0,1), las = 1, cex.axis = 0.8)]
# plotting two additional chlorpromazine replicate populations

chlorpRepsDT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
counter <- 0
lapply(unique(chlorpRepsDT$Replicate), function(p) {
    chlorpRepsDTsmooth <- ksmoother(chlorpRepsDT[Replicate == p])
    chlorpRepsDTsmooth[Replicate == p, points(pos, maxFreq-counter, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
    chlorpRepsDTsmooth[Replicate == p, lines(pos, maxFreq-counter, col = "lightgray", bty = "n", lwd = 0.3)]
    counter <<- counter + 0.1
} )
chlorpRepsDT[, mtext(text=paste0("Chr", as.roman(chrom)), side = 1, line = 1, cex = textSize, xpd = NA)]
dev.off()

# ============================================================================
# Plot the most frequent haplotype at each position at a specific chromosome without y axes 

chrom <- 4
chrDT <- mrgeDT[chr == chrom]
baseDT <- chrDT[Chemical == "Base_population"]
cadRepsDT <- chrDT[Chemical == "cadmium_chloride" & Replicate %in% c("R04", "R10", "R11")]
chlorpRepsDT <- chrDT[Chemical == "chlorpromazine" & Replicate %in% c("R02", "R04", "R10")]

filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/")
textSize <- 1.5
ymax <- 1
ymin <- -0.5

pdf(paste0(filePath, "MBE_cover_submission_v2.pdf"), width = 8, height = 10)
#png(paste0(filePath, "SEE01_all_chems_MIH_hap_devs.png"), width = 8, height = 10, units = 'in', res = 350)
layout(matrix(1:3))
par(mar = c(0,0,0,0), oma = c(0,0,0,0), xpd = NA)
# plotting base population

baseDT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(min(baseDT$maxFreq), max(baseDT$maxFreq)), axes="FALSE", yaxs="i", xaxs="i")]
baseDT[, points(pos, maxFreq, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
#baseDT[, lines(pos, maxFreq, col = "lightgray", bty = "n", lwd = 0.3)]
baseDT[, mtext(text="Ancestral population", side = 3, cex = textSize, xpd = NA, line = -2)]
# plotting three cadmium chloride replicate populations
#par(mar = c(0,0,0,0), oma = c(0,2.5,0,2))
cadRepsDT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(min(cadRepsDT$maxFreq)-0.5, max(cadRepsDT$maxFreq)), axes="FALSE", yaxs="i", xaxs="i")]
counter <- 0
lapply(unique(cadRepsDT$Replicate), function(p) {
    cadRepsDT[Replicate == p, points(pos, maxFreq-counter, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
    #cadRepsDT[Replicate == p, lines(pos, maxFreq-counter, col = "lightgray", bty = "n", lwd = 0.3)]
    counter <<- counter + 0.25
} )
cadRepsDT[, mtext(text="Stressor 1 evolved", side = 3, cex = textSize, xpd = NA, line = -1)]
# plotting three chlorpromazine replicate populations

chlorpRepsDT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(pos),max(pos)), ylim = c(min(chlorpRepsDT$maxFreq)-0.7, max(chlorpRepsDT$maxFreq)), axes="FALSE", yaxs="i", xaxs="i")]
counter <- 0
lapply(unique(chlorpRepsDT$Replicate), function(p) {
    chlorpRepsDT[Replicate == p, points(pos, maxFreq-counter, col = colorCodes, bty = "n", pch = 16, cex = 0.75)]
    #chlorpRepsDT[Replicate == p, lines(pos, maxFreq-counter, col = "lightgray", bty = "n", lwd = 0.3)]
    counter <<- counter + 0.25
} )
chlorpRepsDT[, mtext(text="Stressor 2 evolved", side = 3, cex = textSize, xpd = NA, line = -1)]
chlorpRepsDT[, mtext(text=paste0("Chr", as.roman(chrom)), side = 1, line = -2, cex = textSize, xpd = NA)]
dev.off()
# ============================================================================
# Trouble-shooting