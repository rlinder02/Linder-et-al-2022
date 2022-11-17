#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Finding regions under selection
# Description: Finding regions of significant change between the evolved populations and the base 

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: analysis_step11_calc_abs_snp_difs.R <avg_haplotypes> <helper>", call.=FALSE)
}

haplotypes <- args[1]
helper <- args[2]

# ============================================================================
# Load packages and sourced files

library(tictoc)
library(data.table)
library(R.utils)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

# ============================================================================
# Load data

avgHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "/", analysisType = "hap_freqs_avg_sd_difs", samplePattern = ".*")

# ============================================================================
# Find the top 23 most frequent haplotype combinations across all drugs

findTopHaps <- lapply(avgHapDifsDTs, function(DT) {
    tictoc::tic("Time taken")
    print(DT$chemWeek[1])
    flush.console()
    freqs <- DT[, tstrsplit(avgDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
    founderClasses <- DT[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
    topSixFounders <- unlist(lapply(1:nrow(freqs), function(x) {
        topVals <- freqs[x, order(abs(as.numeric(freqs[x])), decreasing = TRUE)[1:6]]
        topFounders <- founderClasses[x, ..topVals]
        as.character(topFounders)
    } ) )
    toc()
    topSixFounders
} )
allTopHaps <- unlist(findTopHaps)
hapFlavors <- sort(table(allTopHaps),decreasing=TRUE)
topHaps <- paste(names(hapFlavors)[1:23], collapse = ";")
topHapsDT <- data.frame(topHapCombos = topHaps)
write.table(topHapsDT, "topHapCombosAllChems.txt", col.names = T, quote = F)

# ============================================================================
# Trouble-shooting