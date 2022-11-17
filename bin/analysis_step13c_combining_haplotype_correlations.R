#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Combine correlation tables
# Description: Combine the Spearman correlation tables

##############################################################################
# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Usage: analysis_step13c_combining_haplotype_correlations.R <correlations> ", call.=FALSE)
}

correlations <- args[1]

# ============================================================================
# Load packages and sourced files


# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

# ============================================================================
# Load data

directory <- paste0(projectDir, "/")
files <- dir(path = directory, pattern = "^Idx")
setwd(directory)

# ============================================================================
# Run analyses

allCorDFs <- lapply(files, function(read) {
    read.table(read, header = T, sep = "\t")
} )
allCors <- do.call(rbind, allCorDFs)
write.table(allCors, file= "complete_chems_all_reps_spearman_cors.txt", col.names = T, quote = F, sep = "\t")
# ============================================================================
# Trouble-shooting