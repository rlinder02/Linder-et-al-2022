#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Genomewide coverage
# Description: Calculates the mean genomewide coverage per sample, excluding the mitochondria and repetitive regions, and stores this as a table. 

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("Usage: preprocess_step2-coverage_genomewide_calculations_plots.R <snp_cov> <offsets> <repeats> <projectDir>", call.=FALSE)
}

cov <-  args[1]
offsets <- args[2]
repeat_lst <- args[3]
helper <- args[4]

# ============================================================================
# Load packages and sourced files

library(data.table)
library(tictoc)
library(ggpubr)
source(helper)

# ============================================================================
# Set global options

options(digits = 10)
projectDir <- getwd()

# ============================================================================
# Custom functions

# ============================================================================
# Load data

offsets <- fread(offsets, header = T)
g_l <- posoff_chr.bounds(offsets)   
snpCovDT <- fread(cov, header = TRUE)
repeats <- fread(repeat_lst, header = T)
repeats[,cp := paste(C, p, sep = "_")]

# ============================================================================
# Calculate coverages genomewide and write the tables to disk.

plottingBaseCovDT <- covcalc_create.chr.cov.dts(window = 2000, samplePattern = "^BAS02$")
generateBasePlottingCovsDT <- plottingBaseCovDT(as.data.frame(snpCovDT))
plottingCovsDTs <- covcalc_create.chr.cov.dts(window = 2000, samplePattern = "^SEE")
generatePlottingCovsDTs <- plottingCovsDTs(as.data.frame(snpCovDT))
savingPlottingCovsDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "/Tables/Genomewide_plotting_coverages_2kb_window/", analysisType = "plotting_coverage", idCol = "id")
writePlottingBaseCovFile <- savingPlottingCovsDTs(generateBasePlottingCovsDT)
writePlottingCovsFiles <- savingPlottingCovsDTs(generatePlottingCovsDTs)

gwideCovsDT <- covcalc_create.gwide.cov.dts(snpCovDT, samplePattern = "^SEE")
savingGwideCovDT <- writing_write.file(projectRootDir = projectDir, folderPath = "/Tables/Genomewide_avg_coverages/", fileName = "avg_gwide_coverage")
writeGwideCovsFiles <- savingGwideCovDT(gwideCovsDT)

# ============================================================================
# Generate a look-up table of low-coverage samples to filter for downstream analyses; also generate table of coverages for samples using 

covCutoff <- covcalc_filter.low.cov.samples(covCutoff = 5)
lowCovSamples <- covCutoff(gwideCovsDT)
fwrite(lowCovSamples, paste0(projectDir, "/Tables/Genomewide_avg_coverages/low_cov_samples_DT.txt"), sep = "\t")

# ============================================================================
# Generate tables of sample coverages vs base coverages

treatmentDTs <- generatePlottingCovsDTs 
baseDT <- generateBasePlottingCovsDT
baseDTs <- rep(list(baseDT), length(treatmentDTs))
calcNormCovDTs <- Map(covcalc_create.norm.cov.dts, treatmentDTs, baseDTs)
savingNormalizedCovsDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "/Tables/Treatment_vs_base_coverage_tables_2kb/", analysisType = "normalized_coverage", idCol = "id")
writeNormalizedCovsFiles <- savingNormalizedCovsDTs(calcNormCovDTs)

# ============================================================================
# Trouble-shooting
