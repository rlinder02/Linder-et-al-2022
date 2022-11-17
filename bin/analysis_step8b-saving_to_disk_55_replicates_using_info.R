#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Replicates using in paper
# Description: Saving information on the 55 replicates being used throughout this paper

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: analysis_step8b-saving_to_disk_ind_treatments_vs_base_differences.R <transformed_haplotypes> <helper>", call.=FALSE)
}

haplotypes <- args[1]
helper <- args[2]

# ============================================================================
# Load packages and sourced files

library(tictoc)
library(data.table)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()
print(getwd())
flush.console()

# ============================================================================
# Custom functions

# ============================================================================
# Load data


indHapDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "/", analysisType = "haps_sq_diffs", samplePattern = "^SEE12")
allDT <- do.call(rbind, indHapDTs)

# ============================================================================
# Write out to information on the 55 replicates using in this paper

allDT2 <- allDT[, c("chr", "pos", "Chemical", "Replicate", "id", "collapsedFreqs", "collapsedFounders", "heterozygosity", "freqDifs", "freqDifsSq")]
fwrite(allDT2, "Individual_haplotype_differences_55_selected_populations.txt", col.names = T, sep = "\t")

# ============================================================================
# Trouble-shooting