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

inflect <- function(avgDiffDT, colNames, threshold = 1){
    avgDiffDF <- as.data.frame(avgDiffDT, stringsAsFactors = FALSE)
    allCols <- lapply(colNames, function(col) {
        up   <- sapply(1:threshold, function(n) c(avgDiffDF[,col][-(seq(n))], rep(NA, n)))
        down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), avgDiffDF[,col][-seq(length(avgDiffDF[,col]), length(avgDiffDF[,col]) - abs(n) + 1)]))
        a <- cbind(avgDiffDF[,col],up,down)
        minMaxList <- list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
        names(minMaxList) <- c(paste0(col, "_minima"), paste0(col, "_maxima"))
        minMaxList
    } ) 
    allCols
}

add.sig.column <- function(avgDiffDT, filteredMinMaxList, LODcutoff) {
    avgDiffDT$significant <- 0
    minOnly <- unlist(filteredMinMaxList[[1]][1])
    maxOnly <- unlist(filteredMinMaxList[[1]][2])
    minRows <- unlist(minOnly)
    maxRows <- unlist(maxOnly)
    avgDiffDT$significant[minRows] <- 2
    avgDiffDT$significant[maxRows] <- 1
    avgDiffDT[LOD < LODcutoff, significant := 0]
    avgDiffDT
}

lod.support.intervals <- function(LODdrop) {
    function(DT) {
        chromSplit <- split(DT, DT$chr)
        LODintervals <- lapply(chromSplit, function(chromDT) {
            #print(chromDT$chr[1])
            #flush.console()
            chromDT[, row := .I]
            chromDT[, LOD2 := LOD - LODdrop]
            rightLod <- chromDT[chromDT, x.row-i.row, on = .(row > row, LOD <= LOD2), mult="first"]
            leftLod <- chromDT[chromDT, x.row-i.row, on = .(row < row, LOD <= LOD2), mult="last"]
            chromDT$rightBound <- rightLod
            chromDT$leftBound <- leftLod
            if(any(is.na(chromDT$rightBound))) {
                chromDT$rightBound[is.na(chromDT$rightBound)] <- nrow(chromDT) - chromDT$row[is.na(chromDT$rightBound)]}
            if(any(is.na(chromDT$leftBound))) {
                chromDT$leftBound[is.na(chromDT$leftBound)] <- 1 - chromDT$row[is.na(chromDT$leftBound)]}
            chromDT
        } )
        chromDTs <- do.call(rbind, LODintervals)
        chromDTs
    }
}

find.lod.dts <- function(DT, peakGP) {
    DT[, row := .I]
    counter <- 0
    findFlanks <- lapply(peakGP, function(peak) {
        counter <<- counter + 1
        peakRow <- DT[gp == peak, row]
        rowsDown <- DT[gp == peak, rightBound]
        rowsUp <- DT[gp == peak, leftBound]
        lodDT <- DT[row <= (peakRow + rowsDown) & row >= (peakRow + rowsUp)][, Idx := counter]
        lodDT
    } )
    flanksDT <- do.call(rbind, findFlanks)
    flanksDT
}

# ============================================================================
# Load data needed for downstream analyses

founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)
indHapDiffsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_sum_of_squared_diffs_tables/", analysisType = "haps_sq_diffs", samplePattern = "^SEE12B02")
indLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Hap_adjusted_LOD_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^SEE12B02")
indHapDevsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_rep_hap_dev_tables/", analysisType = "avg_hap_devs", samplePattern = "^SEE12B02")
repLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Repeatability_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^18way")

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
# Calculating z-scores for log transformed LOD data and comparing the 3 vs 3 replicates genomewide via Spearmans' rho

popDT <- do.call(rbind, repLODDTs)
popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
lodDT <- dcast(popDT, gp + chr + pos ~ Chemical, value.var = "logLOD")
names(lodDT) <- gsub("18way_", "", names(lodDT))
names(lodDT) <- gsub("_", " ", names(lodDT))
lodDF <- as.data.frame(lodDT)
lodDF[,4:ncol(lodDF)] <- scale(lodDF[,4:ncol(lodDF)])
lodDT <- as.data.table(lodDF)

# ============================================================================
# Plot LOD scores for all pairs of treatments with trend line to look for correlations within each chemical

popDT <- do.call(rbind, repLODDTs)
popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
chemSplit <- split(popDT, popDT$chemWeek)
chemLoop <- lapply(chemSplit, function(chem) {
    print(chem$chemWeek[1])
    flush.console()
    replodDT <- dcast(chem, gp + chr + pos ~ Chemical, value.var = "LOD")
    names(replodDT) <- gsub("18way_", "", names(replodDT))
    names(replodDT) <- gsub("_", " ", names(replodDT))
    replodDF <- as.data.frame(replodDT)
    #replodDF[,4:ncol(replodDF)] <- scale(replodDF[,4:ncol(replodDF)])
    panelPlot <- ggplot(replodDF, aes(replodDF[,4], replodDF[,5])) + xlab("") + ylab("") + facet_grid(names(replodDF)[4] ~ names(replodDF[5]), scales = "free", labeller = labeller(yvar = label_wrap_gen(16))) + geom_hex(bins = 50) + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) 
    panelPlot
} )
savePlot <- do.call(ggarrange, c(chemLoop, list(common.legend = TRUE, legend = "right")))
savingPlot <- annotate_figure(savePlot, bottom = text_grob("LOD", color = "black", size = 9, vjust = -2), left = text_grob("LOD", color = "black", rot = 90, size = 9))
#ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure4A.pdf"), savingPlot, width = 6, height = 5.5, units = "in")
#ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure4A.png"), savingPlot, width = 6, height = 5.5, units = "in", dpi = 350)

# ============================================================================
# Plot super-imposed LOD scores for within drug subsets

myColPal <- add.alpha(mycols2, 0.5)

chroms <- as.character(as.roman(1:16))
popDT <- do.call(rbind, repLODDTs)
popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
chem <- popDT[grepl("cadmium", Chemical)]
replodDT <- dcast(chem, gp + chr + pos ~ Chemical, value.var = "logLOD")
names(replodDT) <- gsub("18way_", "", names(replodDT))
names(replodDT) <- gsub("_", " ", names(replodDT))
chems <- names(replodDT)[4:5]
lodDF <- as.data.frame(replodDT)
lodDF[,4:ncol(lodDF)] <- scale(lodDF[,4:ncol(lodDF)])
replodDT <- as.data.table(lodDF)
replodDT2 <- melt(replodDT, id.vars = c("gp", "chr", "pos"), measure.vars = c("cadmium chloride grp 1", "cadmium chloride grp 2"),variable.name = "group", value.name = "zscore")
replodDT2[, chrBounds := "gray"][(chr %% 2) == 0, chrBounds := "white"]

start <- ch.bounds[seq(1, length(ch.bounds)-2, 2)]
end <- ch.bounds[seq(2, length(ch.bounds)-2, 2)]
chrDT <- data.table(start = start, end = end, chrBounds = rep("gray", 8))

repsvsrepsPlot <- ggplot() + geom_line(data = replodDT2, aes(gp, zscore, colour = group, group = group)) + geom_rect(data = chrDT, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.3) + theme_bw(base_size = 10) + theme(panel.grid = element_blank(), legend.position = "bottom") + scale_colour_manual(values = alpha(c("blue", "red"), 0.3)) + scale_x_continuous(limits = c(0, ch.bounds[17]), breaks = c(g_l[1:16] + offsets[[2]][1:16]/2), labels = as.roman(1:16), expand = c(0, 0)) + coord_cartesian(ylim = c(-4,4)) + xlab("")

#revisedPlot <- reposition_legend(repsvsrepsPlot, 'top right',  offset=0.002)

#ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure4B.pdf"), revisedPlot, width = 6, height = 5.5, units = "in")
#ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure4B.png"), revisedPlot, width = 6, height = 5.5, units = "in", dpi = 350)

fig4 <- ggarrange(savingPlot, repsvsrepsPlot, ncol = 1, labels = "AUTO", heights = c(2,1))

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure6_repeatability_LOD_plots_revised.pdf"), fig4, width = 8, height = 10, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure6_repeatability_LOD_plots_revised.png"), fig4, width = 8, height = 10, units = "in", dpi = 350)

# =====================================================
# Trouble-shooting

popDT <- do.call(rbind, repLODDTs)
popDT <- popDT[id != "SEE12B02YP000"]
popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
chemSplit <- split(popDT, popDT$chemWeek)
chemLoop <- lapply(chemSplit, function(chem) {
    print(chem$chemWeek[1])
    flush.console()
    replodDT <- dcast(chem, gp + chr + pos ~ Chemical, value.var = "LOD")
    names(replodDT) <- gsub("18way_", "", names(replodDT))
    names(replodDT) <- gsub("_", " ", names(replodDT))
    replodDF <- as.data.frame(replodDT)
    #replodDF[,4:ncol(replodDF)] <- scale(replodDF[,4:ncol(replodDF)])
    panelPlot <- ggplot(replodDF, aes(replodDF[,4], replodDF[,5])) + xlab("") + ylab("") + facet_grid(names(replodDF)[4] ~ names(replodDF[5]), scales = "free", labeller = labeller(yvar = label_wrap_gen(16))) + geom_hex(bins = 50) + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) 
    panelPlot
} )
savePlot <- do.call(ggarrange, c(chemLoop, list(common.legend = TRUE, legend = "right")))
savingPlot <- annotate_figure(savePlot, bottom = text_grob("LOD", color = "black", size = 9, vjust = -2), left = text_grob("LOD", color = "black", rot = 90, size = 9))

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure7_repeatability_LOD_plots_v2.tiff"), device = "tiff", savingPlot, width = 8, height = 6, units = "in", dpi = 350)
