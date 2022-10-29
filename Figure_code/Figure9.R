##############################################################################

# Project:     DXQTL02
# Author:      Robert Linder
# Date:        2019
# Title:       title
# Description: A brief description

##############################################################################

# ============================================================================
# Load packages and sourced files
# Sourced files are kept in the default working directory	
library(rtracklayer)
library(Gviz)
library(R.utils)
library(GenomicRanges)
library(tictoc)
library(data.table)
library(tidyverse)
library(scales)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggbeeswarm)
library(ggExtra)
library(DEGreport)
library(DescTools)
library(GISTools)
library(pcaPP)
library(grid)
library(gridExtra)
library(cowplot)
library(gtable)
library(NbClust)
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

makePairs <- function(data) 
{
    grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
    #grid <- subset(grid, x != y)
    all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
        xcol <- grid[i, "x"]
        ycol <- grid[i, "y"]
        data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
                   x = data[, xcol], y = data[, ycol], data)
    }))
    all$xvar <- factor(all$xvar, levels = names(data))
    all$yvar <- factor(all$yvar, levels = names(data))
    densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
        data.frame(xvar = names(data)[i], yvar = names(data)[i], x = data[, i])
    }))
    list(all=all, densities=densities)
}

makeUniquePairs <- function(data) 
{
    grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
    grid <- subset(grid, x != y)
    gridDT <- as.data.table(grid)
    dt2 <- unique(gridDT, by=c("x", "y"))[x > y, c("x", "y") := .(y, x)]
    res2 <- unique(dt2, by=c("x", "y"))
    res3 <- as.data.frame(res2)
    all <- do.call("rbind", lapply(1:nrow(res3), function(i) {
        xcol <- res3[i, "x"]
        ycol <- res3[i, "y"]
        data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
                   x = data[, xcol], y = data[, ycol], data)
    }))
    all$xvar <- factor(all$xvar, levels = names(data))
    all$yvar <- factor(all$yvar, levels = names(data))
    densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
        data.frame(xvar = names(data)[i], yvar = names(data)[i], x = data[, i])
    }))
    list(all=all, densities=densities)
}

allUniquePairs <- function(data) 
{
    grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
    grid <- subset(grid, x != y)
    gridDT <- as.data.table(grid)
    dt2 <- unique(gridDT, by=c("x", "y"))[x > y, c("x", "y") := .(y, x)]
    res2 <- unique(dt2, by=c("x", "y"))
    colCombos <- data.frame(apply(res2, 1, function(x) names(data)[x]))
    colCombos
}

# ============================================================================
# Load data needed for downstream analyses

privateSNPDT <- fread('Private_snps_AB3B6_lookup_table.txt', header = T)
privateSNPDT$unique_founder[grep("_", privateSNPDT$unique_founder)] <- "AB3"
avgSnpDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Avg_snps_abs_diffs_dfs/", analysisType = "snps_abs_diffs", samplePattern = ".*")
avgSnpDifsDTs2 <- lapply(avgSnpDifsDTs, function(x) {
    x$founder <- NA
    x$founder[x$gp %in% privateSNPDT$gp] <- privateSNPDT$unique_founder
    names(x)[5] <- "control.difs"
    x
})

indLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Hap_adjusted_LOD_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^SEE12B02")
repLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Repeatability_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^18way")
avgHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Averaged_and_sd_tx_tables/", analysisType = "hap_freqs_avg_sd_difs", samplePattern = ".*")

treatKeyDT <- fread("treatment_key.txt", header = T)
topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))
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
mypal <- RColorBrewer::brewer.pal(4, "Dark2")
myColPal <- add.alpha(mycols, 0.5)
treatmentAbbrv <- fread("chemAbbrev.txt")

# ============================================================================
# Plot log-transformed LOD scores for all pairs of treatments with trend line to look for correlations b/w treatments

popDT <- do.call(rbind, indLODDTs)
popDT[, Chemical := treatmentAbbrv$abbr[treatmentAbbrv$chemical3 == Chemical], by = Chemical]
lodDT <- dcast(popDT, gp + chr + pos ~ Chemical, value.var = "LOD")
names(lodDT) <- gsub("18way_", "", names(lodDT))
names(lodDT) <- gsub("_", " ", names(lodDT))
#allCor <- cor(lodDT[, 4:ncol(lodDT)], method = "kendall")
lodDF <- as.data.frame(lodDT)
gg1 = makePairs(lodDF[,-c(1:3)])
mega_iris = data.frame(gg1$all)
megaIrisDT <- as.data.table(mega_iris)
maxVals <- megaIrisDT[, .(x = max(x), y = max(y)), by = c("xvar", "yvar")][, c("x", "y") := .(1.2*x, 1.2*y)]
fileName <- "Figure9_correlated_LOD_scores_plot.pdf"
        
g <- ggplot(mega_iris, aes_string(x = "x", y = "y")) + xlab("") + ylab("") + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "center", label.y.npc = "top", size = 3, hjust = 0.5) + geom_hex(bins = 50) + facet_grid(xvar ~ yvar, scales = "free", labeller = labeller(yvar = label_wrap_gen(8), xvar = label_wrap_gen(8))) + scale_fill_continuous(type = "viridis", ) + geom_smooth(span = 0.3) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1), aspect.ratio = 1) + scale_x_continuous(limits = c(0, NA), expand = c(0,0)) + scale_y_continuous(limits = c(0, NA), expand = c(0,0)) + xlab("-log10(p)") + ylab("-log10(p)")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/", fileName), g, width = 8, height = 10, units = "in")
system2(command = "pdfcrop", args = c(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure9_correlated_LOD_scores_plot.pdf"), paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure9_correlated_LOD_scores_plot_cropped.pdf")))

# ============================================================================
# Plot LOD scores for all pairs of treatments with trend line to look for correlations b/w treatments WITHOUT chr X; doesn't really change correlations...

popDT <- do.call(rbind, indLODDTs)
popDT2 <- popDT[chr != 10]
lodDT <- dcast(popDT2, gp + chr + pos ~ Chemical, value.var = "LOD")
names(lodDT) <- gsub("18way_", "", names(lodDT))
names(lodDT) <- gsub("_", " ", names(lodDT))
#allCor <- cor(lodDT[, 4:ncol(lodDT)], method = "kendall")
lodDF <- as.data.frame(lodDT)
gg1 = makePairs(lodDF[,-c(1:3)])
mega_iris = data.frame(gg1$all)
megaIrisDT <- as.data.table(mega_iris)
maxVals <- megaIrisDT[, .(x = max(x), y = max(y)), by = c("xvar", "yvar")][, c("x", "y") := .(1.2*x, 1.2*y)]
fileName <- "Figure6v2.pdf"

g <- ggplot(mega_iris, aes_string(x = "x", y = "y")) + xlab("") + ylab("") + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "center", label.y.npc = "top", size = 3, hjust = 0.5) + geom_hex(bins = 50) + facet_grid(xvar ~ yvar, scales = "free", labeller = labeller(yvar = label_wrap_gen(8), xvar = label_wrap_gen(8))) + scale_fill_continuous(type = "viridis", ) + geom_smooth(span = 0.3) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) + scale_x_continuous(limits = c(0, NA), expand = c(0,0)) + scale_y_continuous(limits = c(0, NA), expand = c(0,0))
savingPlot <- annotate_figure(g, bottom = text_grob("-log10(p)", color = "black", size = 14), left = text_grob("-log10(p)", color = "black", rot = 90, size = 14))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/", fileName), savingPlot, width = 8, height = 10, units = "in")



# ============================================================================
# Trouble-shooting
