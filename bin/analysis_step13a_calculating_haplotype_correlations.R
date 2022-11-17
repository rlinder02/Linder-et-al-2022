#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Haplotype frequency correlations
# Description: Preparing the files needed for calculating the Spearman correlation genome-wide of haplotype frequencies across all samples not eliminated from further analyses based on heterozygosity

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: analysis_step13a_calculating_haplotype_correlations.R <transformed_haplotypes> <helper>", call.=FALSE)
}

hap_diffs <-  args[1]
metadata <- args[2]
helper <- args[3]

# ============================================================================
# Load packages and sourced files

library(tictoc)
library(data.table)
library(tidyverse)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

# ============================================================================
# Load data 

indHapDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "/Ind_tx_ind_haps_het_diffs_tables/", analysisType = "haps_het_diffs", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]_")
allReps <- fread(metadata, header = T) 

# ============================================================================
# Checking to ensure have same vector of haplotypes at a given position between chemical treatments (there are differences between chemicals, but not within, so test is valid) for ALL chemicals

popDTs <- do.call(rbind, indHapDTs)
collHaps <- popDTs[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
collHapsDT <- cbind(popDTs[, "gp"], collHaps)
gpSplit <- split(collHapsDT, collHapsDT$gp)
gpLooper <- lapply(gpSplit, function(gpos) {
    checkEquality <- apply(gpos, 2, function(x) length(unique(x)) == 1)
    if(any(checkEquality == FALSE)) {return(gpos$gp[1])}
} )
anygp <- do.call(rbind, gpLooper)
uneqGPsAll <- sort(unique(anygp))

# ============================================================================
# Calculating pairwise correlations amongst ALL replicates across ALL treatments not eliminated based on heterozygosity

popDTs <- do.call(rbind, indHapDTs)
popDTs <- popDTs[!gp %in% uneqGPsAll]
sexualReps <- allReps[Type == "Outbred_sexual"]
popDTs <- popDTs[id %in% sexualReps$id]
PopDTcpy <- copy(popDTs)
allGPsDT <- data.table(gp = unique(PopDTcpy$gp))
allGPsDT[, Idx := c(0, rep(1:(nrow(allGPsDT)-1)%/%1000))] ### Can split into 11 groups to parallelize on the hpcc
PopDTcpyMrge <- merge(PopDTcpy, allGPsDT, by = "gp")
popDTs <- merge(popDTs, allGPsDT, by = "gp")
popDT2 <- popDTs[, ID := paste0(chemWeek, "-",  gp, ".", Idx, "-", Replicate)][, c("Chemical", "chemWeek", "Replicate", "cp", "chr", "pos", "id", "gp", "Week") := NULL]
freqs <- popDT2[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
popDT2 <- cbind(freqs, popDT2[, "ID"])
popDF <- as.data.frame(popDT2, stringsAsFactors = FALSE)
popDF2 <- t(popDF)
popDF2 <- as.data.frame(popDF2, stringsAsFactors = FALSE)
names(popDF2) <- unique(popDF$ID)
popDF2 <- popDF2[-nrow(popDF2),]
popDF3 <- as.data.frame(apply(popDF2, 2, as.numeric), stringsAsFactors = FALSE)
popDT3 <- as.data.table(popDF3)
chemWeekPos <- unique(paste0(PopDTcpyMrge$chemWeek, "-", PopDTcpyMrge$gp, ".", PopDTcpyMrge$Idx))
idces <- unique(PopDTcpyMrge$Idx)
idxDT <- as.data.table(idces)
fwrite(idxDT, "SEE01_gposIdx.txt", col.names = F)
fwrite(popDF3, "Reformatted_SEE01_hap_freqs_v3.txt")
fwrite(as.data.table(chemWeekPos), "Unique_chemPos_v2.txt", col.names = FALSE)

# ============================================================================
# Trouble-shooting