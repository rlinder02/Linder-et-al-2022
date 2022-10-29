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


# ============================================================================
# Calculating average per-site correlation

# directory <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/")
# files <- dir(path = directory, pattern = "_spearman_cors.txt$")
# setwd(directory)
# #allCorDFs <- lapply(files, function(read) {
# #read.table(read, header = T, sep = "\t")
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

# ============================================================================
# Plotting avg heterozygosity vs avg correlation per treatment 

corJob <- read.table("sexual_reps_6_plus_spearman_avgs.txt", header = T, sep = "\t")
corJob <- corJob[-grep("fluconazole|nic", rownames(corJob)),]
corJob <- corJob[-grep("cadmium_chloride.*-R0[1-3]|YPD.*-R10|glacial_acetic_acid.*-R08|chlorpromazine.*-R12|chlorpromazine.*-R13", rownames(corJob)),]
corJob <- corJob[, -grep("fluconazole|nic|cadmium_chloride.R0[1-3]|YPD.R10|glacial_acetic_acid.R08|chlorpromazine.R12|chlorpromazine.R13", names(corJob))]
chemReps <- unique(names(corJob))
chemNames <- unique(gsub(".R[0-9][0-9]", "", chemReps))

hapDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"), header = T)
hapDT$treatment <- gsub("^18way_|_([0-1]?[0-9])$", "", hapDT$chemical)
avgHets <- hapDT[, mean(het), by = chemical]
setnames(avgHets, c("chemWeek", "avgHet"))
avgHets$id <- unlist(lapply(avgHets$chemWeek, function(x) {
    splitting <- strsplit(x,"_")[[1]]
    if(length(splitting) == 3) {return(splitting[2])} else{
        return(paste(splitting[2:3], collapse = "_"))}
} ) )

avgChemCorDFs <- lapply(unique(chemNames), function(x) {
    chemDF <- corJob[grepl(x, rownames(corJob)), grepl(x, names(corJob))]
    avgCor <- mean(chemDF[upper.tri(chemDF)])
    avgCorDF <- data.frame(chemical = x, avgCor = avgCor, stringsAsFactors = FALSE)
    avgCorDF
} )


# avgChemCorDFs <- lapply(corJob$chemLoop, function(x) {
#     x <- as.data.frame(x)
#     rownames(x) <- x$id
#     names(x)[2:length(names(x))] <- x$id
#     chemName <- gsub("-.*", "", names(x)[2])
#     x$id <- NULL
#     avgCor <- mean(x[upper.tri(x)])
#     avgCorDF <- data.frame(chemical = chemName, avgCor = avgCor, stringsAsFactors = FALSE)
#     avgCorDF
# } )

avgCorsDF <- do.call(rbind, avgChemCorDFs)
avgCorDT <- as.data.table(avgCorsDF)
avgHets[, chemical := gsub("18way_|_12", "", chemWeek)][, chemWeek := NULL]
avgCorHetDF <- merge(avgHets, avgCorDT, by = "chemical")
avgCorHetDF[chemical %like% "YPD", c("chemical", "id") := .('ypd', 'ypd')]
avgCorHetDF <- avgCorHetDF[order(avgCorHetDF$chemical),]
avgCorHetDF <- na.omit(avgCorHetDF)
avgCorHetDF$chemical[7] <- "YPD"

chemsUsing <- gsub("_", " ", avgCorHetDF$chemical)
chemsUsing[chemsUsing == "YPD"] <- "ypd"
avgCorHetDF$chemical <- gsub("_", " ", avgCorHetDF$chemical)


gplot <- ggplot(data=avgCorHetDF, aes(avgHet, avgCor)) + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(aes(label = ..rr.label..)) + geom_point(aes(colour = chemical), size = 2) + xlab("average within-treatment heterozygosity") + ylab("average within-treatment correlation") + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + scale_colour_manual(values = colorPal[names(colorPal) %in% avgCorHetDF$chemical])

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS14.pdf"), gplot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS14.png"), gplot, width = 8, height = 9, units = "in", dpi = 350)

# ============================================================================
# Trouble-shooting