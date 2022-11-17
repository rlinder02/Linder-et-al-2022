#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Reformatting haplotype data
# Description: Reformat the haplotype data so it is easier to work with.

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("Usage: preprocess_step3-reformatting_haplotype_call_data.R <haplotype_frequencies> <founder_names> <offsets> <helper>", call.=FALSE)
}

haplotypes <-  args[1]
founders <- args[2]
offsets <- args[3]
helper <- args[4]

# ============================================================================
# Load packages and sourced files

library(data.table)
library(R.utils)
source(helper)

# ============================================================================
# Set global options

# ============================================================================
# Custom functions

# ============================================================================
# Load data

offsets <- fread(offsets, header = T)
g_l <- posoff_chr.bounds(offsets)
hapDT <- fread(haplotypes,header=TRUE)
founderNames <- fread(founders, header = F)
founderNames <- founderNames[,V1]

# ============================================================================
# Reformat and save reformatted haplotype frequency table

colsKeeping <- c("poolroot", "chr", "pos", "cutree", "founderfreqs")
sample.pattern <- hapfilereform_reformat.haplotype.file(samplePattern = "^SEE", basePop = "^BAS02$|^DIP02", colsKeeping = colsKeeping, idCol = "poolroot")
reformatHapDT <- sample.pattern(hapDT)
checkHapDT <- hapfilereform_check.haplotype.file(reformatHapDT, idCol = "poolroot")
fwrite(reformatHapDT,"Jan2022.allhaps.restructured.txt", quote = F, sep = "\t", col.names = T, row.names = F)
gzip("Jan2022.allhaps.restructured.txt", destname = "Jan2022.allhaps.restructured.txt.gz")

# ============================================================================
# Trouble-shooting
