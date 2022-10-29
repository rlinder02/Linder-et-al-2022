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
library(lemon)
library(abind)
library(job)
library(ggdendro)
library(cowplot)
library(magick)
library(ggh4x)
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
indHapDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_het_diffs_tables/", analysisType = "haps_het_diffs", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]_")
indHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_het_diffs_tables/", analysisType = "haps_het_diffs", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]_")
allBind <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_All_Reps_Hets.txt"))
topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))
hetDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"), header = T)
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
mycols2 <- c("blue", "red")
myColPal <- add.alpha(mycols2, 0.5)

# ============================================================================
# Data wrangle frequency change of all haplotypes (for each replicate), with each replicate a different shape and the haplotypes colored consistently with previous plots.

popDT <- do.call(rbind, indHapDTs)
popDT2 <- popDT[id %in% c("SEE12B02CD600R04", "SEE12B02GA120R15", "SEE12B02CA018R10")]
chemSplit <- split(popDT2, popDT2$id)

chemLooper <- lapply(chemSplit, function(chem) {
    freqs <- chem[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
    collHaps <- chem[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
    allHaps <- cbind(chem$chr, chem$pos, chem$gp, chem$Chemical, chem$Replicate, chem$id, collHaps, freqs)
    freqsDF <- as.data.frame(freqs)
    newFreqs <- lapply(freqsDF, function(y) gsub("^-\\d.*|^\\d.*", 1, y))
    newFreqsDF <- as.data.frame(do.call(cbind, newFreqs))
    newFreqsDF2 <- sapply(newFreqsDF, as.numeric)
    allHaps$reps <- rowSums(newFreqsDF2, na.rm = T)
    positions <- rep(unlist(allHaps[,2]), allHaps$reps) 
    gpositions <- rep(unlist(allHaps[,3]), allHaps$reps)
    replicates <- sort(rep(unlist(allHaps[,5]), allHaps$reps))
    allfreqs <- as.numeric(unlist(t(freqs)))
    allfreqs2 <- na.omit(allfreqs)
    hapVec <- as.character(unlist(t(collHaps)))
    hapVec2 <- na.omit(hapVec)
    idxDT <- data.table(chr = allHaps$V1[1], pos = positions, gp = gpositions, Chemical = allHaps$V4[1], Replicate = replicates, id = allHaps$V6[1], haps = hapVec2, freqs = allfreqs2)
    idxDT[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
    hapColors <- match(idxDT$haps, topHaps)
    hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
    Colors <- mycols[hapColors]
    idxDT[, "color.codes" := .(Colors)]
    chem <- idxDT[, hapGrps := haps][!hapGrps %in% topHaps, hapGrps := "other"][order(haps, gp)][, c("difCol") := .(c(1000, diff(gp))), by = "haps"][, row := .I]
    chem[, grpCol := paste0(haps, "_", difCol)]
    chem2 <- chem[, grp := cumsum(c(TRUE, diff(gp)!=1000)), by = haps][order(row)]
    chem2[, grp := cumsum(c(TRUE, diff(gp)!=1000)), by = haps]
    chem2[, grp2 := rleid(haps)]
    chem2[, grpCor := c(0, as.numeric((diff(grp) != 0)))]
    chem2[, grpCor2 := c(0, as.numeric((diff(grp2) != 0)))]
    chem2[, grpSum := rowSums(chem2[, c("grpCor", "grpCor2")])]
    regrp <- chem2[grpSum == 2, grpSum := 1]
    grpThrough <- rle(regrp$grpSum)
    grpPass <- data.table(lengths = grpThrough$lengths, values = grpThrough$values)[, row := .I]
    grpCnt <- grpPass[lengths > 1 & values == 1]
    if(nrow(grpCnt) > 0) {
        newCnt <- lapply(1:nrow(grpCnt), function(x){
            newVals <- 1:grpCnt$lengths[x]
            newRows <- seq(grpCnt$row[x], grpCnt$row[x] + length(newVals)/100, length.out = length(newVals)+1)
            newRows <- newRows[-length(newRows)]
            newDT <- data.table(lengths = 1, values = newVals, row = newRows)
            newDT
        } )
        newCntDT <-  do.call(rbind, newCnt)
        allCntDT <- rbind(grpPass, newCntDT)
    } else {
        allCntDT <- grpPass
    }
    allCnts <- allCntDT[!duplicated(row, fromLast = TRUE)][order(row)]
    allCnts[, grp := rleid(values)]
    allCnts2 <- allCnts[values == 0, grp := (grp-1)]
    finalGrps <- rep(allCnts2$grp, allCnts2$lengths)
    chem2[, finalGrp := finalGrps][, Class := "class"]
    chem2
} )

# ============================================================================
# Plot the frequency of all haplotypes (for each evolved class), with the haplotypes colored consistently with previous plots.

allIdces <- do.call(rbind, chemLooper)
allIdces[Chemical == "caffeine", class := "Aneuploid haploid"][Chemical == "cadmium chloride", class := "Outbred sexual"][Chemical == "glacial acetic acid", class := "Clonal diploid"]
addHet <- merge(allIdces, popDT2[, c("heterozygosity", "gp", "id")], by = c("gp", "id"))
addHet[, chrBounds := "gray"][(chr %% 2) == 0, chrBounds := "white"]
avgHet <- addHet[, .(het = mean(heterozygosity)), by = "class"]
avgHet[, hetLabel := sprintf("het == %.2f", het)]
chroms <- as.character(as.roman(1:16))

start <- ch.bounds[seq(1, length(ch.bounds)-2, 2)]
end <- ch.bounds[seq(2, length(ch.bounds)-2, 2)]
chrDT <- data.table(start = start, end = end, chrBounds = rep("gray", 8))

hapPlots <- ggplot() + geom_point(data = addHet, aes(gp, freqs, colour = hapGrps, group = finalGrp), size = 0.2) + geom_rect(data = chrDT, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.3) + scale_x_continuous(limits = c(0, ch.bounds[17]), breaks = c(g_l[1:16] + offsets[[2]][1:16]/2), labels = as.roman(1:16), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), expand = c(0, 0)) + facet_wrap(~class, ncol = 1) + xlab("") + ylab("haplotype frequency") + scale_colour_manual(values=setNames(addHet$color.codes, addHet$hapGrps)) + theme_bw(base_size = 10) + theme(panel.grid = element_blank(), legend.position = "bottom") + geom_text(data = avgHet, aes(x = Inf, y = Inf, hjust = 1.25, vjust = 1.25, label = hetLabel), parse = TRUE, inherit.aes = FALSE) 
legendPlot <- get_legend(hapPlots + guides(color = guide_legend(nrow = 3, override.aes = list(size = 2), title = NULL)) + theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin()))
hapPlots2 <- hapPlots + theme(legend.position="none")

# ============================================================================
# Plotting histograms of per-site heterozygosity per treatment per replicate with average per-site difference from base

popDT <- do.call(rbind, indHapDifsDTs)
popDT$treatment <- sub("_", "-", popDT$Chemical)
popDT$treatment <- gsub("^(.*)-", "", popDT$treatment)
popDT2 <- popDT[id %in% c("SEE12B02CD600R04", "SEE12B02GA120R15", "SEE12B02CA018R10")]
popDT2[, class := "class"]
popDT2[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
popDT2[Chemical == "caffeine", class := "Aneuploid haploid"][Chemical == "cadmium chloride", class := "Outbred sexual"][Chemical == "glacial acetic acid", class := "Clonal diploid"]

avgHet <- popDT2[, .(het = mean(heterozygosity)), by = class]
avgHet[, hetLabel := sprintf("het == %.2f", het)]

hetPlot <- ggplot() + geom_histogram(data=popDT2, aes(x = heterozygosity), binwidth = 0.01) + geom_text(data = avgHet, aes(x = Inf, y = Inf, hjust = 1.25, vjust = 1.25, label = hetLabel), parse = TRUE, inherit.aes = FALSE) + facet_wrap(~class, ncol = 1) + theme_bw(base_size = 10) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = -45)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 2000))

# ============================================================================
# Putting both together

savePlot <- plot_grid(hapPlots2, hetPlot, align = 'vh',labels = c("A", "B"), ncol = 2, hjust = -1, rel_widths = c(1,0.25))
legendAdjust <- plot_grid(legendPlot, NULL, nrow = 1, rel_widths = c(1, 0.25))
addLegend <- plot_grid(savePlot, legendAdjust, ncol = 1, rel_heights = c(1, .1))

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure2_population_types.pdf"), addLegend, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure2_population_types.png"), addLegend, width = 8, height = 9, units = "in", dpi = 350) 

### Below is for resizing images, keeping the current aspect ratio- unfortunately, images are blurry when shrunk; better to do it manually

imageLoad <- image_read(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure2_population_types.png"))

resizeImage <- image_resize(imageLoad, "502x446")
image_write(resizeImage, path = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure2_population_types_resized2.png"), format = "png")

# ============================================================================
# Trouble-shooting
+ force_panelsizes(rows = unit(3, "in"), cols = unit(4, "in"))
