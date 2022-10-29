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

rep.list <- function(object, repObject) {
    rep(list(object), length(repObject))
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

# ============================================================================
# Plot a pairwise correlation matrix of all replicates for each sample using Spearman correlation

# directory <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/")
# files <- dir(path = directory, pattern = "_spearman_cors.txt$")
# setwd(directory)
# #allCorDFs <- lapply(files, function(read) {
#     #read.table(read, header = T, sep = "\t")
# #} )
# allCorDF <- read.table("all_chems_all_reps_spearman_cors.txt", header = T, sep = "\t")
# setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))
# 
# allCorDF <- allCorDF[-grep("fluconazole|nic", rownames(allCorDF)),]
# allCorDF <- allCorDF[-grep("cadmium_chloride_12.*-R0[1-3]|YPD_12.*-R10|glacial_acetic_acid_12.*-R08|chlorpromazine_12.*-R12|chlorpromazine_12.*-R13", rownames(allCorDF)),]
# allCorDF <- allCorDF[, -grep("fluconazole|nic|cadmium_chloride_12.R0[1-3]|YPD_12.R10|glacial_acetic_acid_12.R08|chlorpromazine_12.R12|chlorpromazine_12.R13", names(allCorDF))]
# 
# job::job(corJob = {
#     chemReps <- unique(names(allCorDF))
#     chemNames <- gsub("X18way_|_12.*", "", chemReps)
#     names(chemReps) <- chemNames
#     chemSplit <- split(chemReps, names(chemReps))
#     chemLoop <- lapply(chemSplit, function(chem) {
#         avgLoop <- lapply(chem, function(x) {
#             tic()
#             corrGrp <- allCorDF[,grep(x, names(allCorDF))]
#             corrGrp <- as.data.frame(corrGrp)
#             names(corrGrp) <- x
#             rownames(corrGrp) <- rownames(allCorDF)
#             corrGrp$id <- gsub(".*18way_", "", rownames(allCorDF))
#             corrGrp$id <- gsub("_12.*-", "-", corrGrp$id)
#             #corrGrp$gp <- as.numeric(gsub(".*_12-|-R.*", "", rownames(allCorDF)))
#             corrGrpDT <- as.data.table(corrGrp)
#             chemName <- gsub("X18way_|_12.*", "", names(corrGrp)[1])
#             corrGrpDT <- corrGrpDT[grepl(chemName, id)]
#             corrAvgs <- corrGrpDT[, mean(get(x), na.rm = T), by = id]
#             setnames(corrAvgs, c("id", x))
#             toc()
#             corrAvgs
#         } )
#         chemAvgs <- Reduce(function(x,y) merge(x = x, y = y, by = "id"), avgLoop)
#         chemAvgs
#     } )
# }, import = c(allCorDF), packages = c("data.table", "tictoc") )

###  Previously generated the necessary table - use this below

corJob <- read.table("sexual_reps_6_plus_spearman_avgs.txt", header = T, sep = "\t")
corJob <- corJob[-grep("fluconazole|nic", rownames(corJob)),]
corJob <- corJob[-grep("cadmium_chloride.*-R0[1-3]|YPD.*-R10|glacial_acetic_acid.*-R08|chlorpromazine.*-R12|chlorpromazine.*-R13", rownames(corJob)),]
corJob <- corJob[, -grep("fluconazole|nic|cadmium_chloride.R0[1-3]|YPD.R10|glacial_acetic_acid.R08|chlorpromazine.R12|chlorpromazine.R13", names(corJob))]
chemReps <- unique(names(corJob))
chemNames <- unique(gsub(".R[0-9][0-9]", "", chemReps))

plotLooper <- lapply(unique(chemNames), function(x) {
    chemDF <- corJob[grepl(x, rownames(corJob)), grepl(x, names(corJob))]
    names(chemDF) <- gsub(".*\\.", "", names(chemDF))
    rownames(chemDF) <- gsub(".*-", "", names(chemDF))
    data <- as.matrix(chemDF)
    df <- reshape2::melt(data)
    gplot <- ggplot(data = df) + geom_tile(aes(x = factor(Var1), y = factor(Var2), fill = value)) + theme_bw(base_size = 8) + scale_fill_continuous(limits=c(0, 1), low = "white", high = "darkred") +  theme(text = element_text(size=10)) + theme(axis.title=element_blank()) + theme(legend.title=element_blank()) + ggtitle(x) + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank()) + theme(axis.text.x = element_text(angle = 45)) + coord_equal()
    gplot
} )
arrangePlot <- do.call(ggarrange, c(plotLooper, list(common.legend = TRUE, legend = "right")))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/FigureS12.pdf"), arrangePlot, width = 8.5, height = 8.5, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/FigureS12.png"), arrangePlot, width = 8.5, height = 8.5, units = "in")


# plotLooper <- lapply(corJob$chemLoop, function(x) {
#     CC3 <- as.data.frame(x)
#     rownames(CC3) <- x$id
#     CC4 <- CC3[, -grep("^id", names(CC3))]
#     names(CC4) <- gsub("X18way_", "", names(CC4))
#     names(CC4) <- gsub("_12.", "-", names(CC4))
#     chemTitle <- gsub("-.*", "", names(CC4))[1]
#     names(CC4) <- gsub(".*-", "", names(CC4))
#     rownames(CC4) <- gsub(".*-", "", rownames(CC4))
#     data <- as.matrix(CC4)
#     df <- reshape2::melt(data)
#     gplot <- ggplot(data = df) + geom_tile(aes(x = factor(Var1), y = factor(Var2), fill = value)) + theme_bw(base_size = 8) + scale_fill_continuous(limits=c(0, 1), low = "white", high = "darkred") +  theme(text = element_text(size=10)) + theme(axis.title=element_blank()) + theme(legend.title=element_blank()) + ggtitle(chemTitle) + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank()) + theme(axis.text.x = element_text(angle = 45)) + coord_equal()
#     gplot
# } )
# arrangePlot <- do.call(ggarrange, c(plotLooper, list(common.legend = TRUE, legend = "right")))
#savePlot <- annotate_figure(arrangePlot, top = text_grob("SEE01_wk12_replicate_correlations", color = "black", size = 10, hjust = 0.5, vjust = -0.05))

# ============================================================================
# Trouble-shooting