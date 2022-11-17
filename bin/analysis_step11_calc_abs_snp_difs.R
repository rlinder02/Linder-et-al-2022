#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       SNP frequency differences
# Description: Calculate the absolute SNP frequency differences between evolved populations and the base population

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: analysis_step11_calc_abs_snp_difs.R <snp_frequencies> <founders> <helper>", call.=FALSE)
}

snps <- args[1]
founders <- args[2]
helper <- args[3]

# ============================================================================
# Load packages and sourced files

library(tictoc)
library(data.table)
library(preprocessCore)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

# ============================================================================
# Load data 

founderNames <- fread(founders, header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)
snpDT <- fread(snps, header=TRUE)

# ============================================================================
# Calculate the absolute differences between the avg snp frequencies of evolved populations and bas02

snpDF <- as.data.frame(snpDT)
samples <- names(snpDF)[grepl("^SEE.*[0-9][0-9][0-9]$", names(snpDF)) | grepl("BAS02", names(snpDF))]
samples <- samples[-grep("SEE12000YP000", samples)]
treatmentSamples <- samples[-grep("BAS02", samples)]
baseSample <- samples[grep("^BAS02$", samples)]
mathFUN <- "abs"
calculateAbsSnpDifsDTs <- Map(calcstats_snp.dif.transforms, treatmentSamples, rep(list(baseSample), length(treatmentSamples)), rep(list(mathFUN), length(treatmentSamples)), rep(list(snpDT), length(treatmentSamples)))
savingSnpDifsDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "/Tables/Avg_snps_abs_diffs_dfs/", analysisType = "snps_abs_diffs", idCol = "id")
writeSnpDifsDTs <- savingSnpDifsDTs(calculateAbsSnpDifsDTs)

# ============================================================================
# Trouble-shooting