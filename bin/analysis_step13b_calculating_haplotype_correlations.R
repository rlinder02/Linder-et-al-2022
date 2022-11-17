#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Haplotype frequency correlations
# Description: Calculating the Spearman correlation genome-wide of haplotype frequencies across all samples not eliminated from further analyses based on heterozygosity

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: analysis_step13b_calculating_haplotype_correlations.R <haplotype_frequencies> <gpositions> <chemPos>", call.=FALSE)
}

gpositions <-  args[1]
hapDT <- args[2]
chemPos <- args[3]

# ============================================================================
# Load packages and sourced files

library(tictoc)
library(data.table)
library(pcaPP)

# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

# ============================================================================
# Load data 

hapFreqs <- fread(hapDT, header = T)
uniquePos <- fread(chemPos, header = F)
nSexualReps <- 105
gpositions <- gsub("\\[", "", gpositions)
gpositions <- gsub("\\]", "", gpositions)

hapFreqsFilt <- hapFreqs[, grep(paste0("\\.", gpositions), colnames(hapFreqs)), with = FALSE]
newCols <- gsub("-R[0-9][0-9]$", "", colnames(hapFreqsFilt))
chemPOS <- data.table(gpIdx = uniquePos[V1 %like% paste0("\\.", gpositions)])
chemPOS[, gp := gsub(".*-", "", gpIdx.V1)][, gp := gsub(paste0("\\.", gpositions), "", gp)]
gpSplit <- split(chemPOS, chemPOS$gp)

tic("Total time")
print(gpositions)
flush.console()
listofCors <- lapply(gpSplit, function(x) {
    if(length(x) > 2) {print(x)}
    findPos <- which(newCols %in% x$gpIdx.V1)
    corrDT <- hapFreqsFilt[, c(findPos), with = FALSE]
    names(corrDT) <- gsub(paste0("\\.", gpositions), "", names(corrDT))
    corrDT <- na.omit(corrDT)
    cors <- cor(corrDT, use = "everything", method = "spearman")
    #cors <- cor.fk(corrDF)
    corsDF <- as.data.frame(cors, stringsAsFactors = FALSE)
    names(corsDF) <- unlist(lapply(names(corsDF), function(x) {paste(strsplit(x, "-")[[1]][c(1,3)], collapse = "-")}))
    if(ncol(corsDF) != nSexualReps) {
        print(c(x$gp[1], ncol(corsDF)))
        flush.console()
        stop("Wrong number of columns")}
    corsDF
})
allCors <- as.data.frame(do.call(rbind, listofCors), stringsAsFactors = FALSE)
toc()
write.table(allCors, file=paste0(projectDir,"/Tables/Correlation_tables/Idx_", gpositions, "_spearman_cors.txt"), col.names = T, quote = F, sep = "\t")