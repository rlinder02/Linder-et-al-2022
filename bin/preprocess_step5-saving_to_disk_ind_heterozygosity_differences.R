#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Calculate haplotype and heterozygosity differences
# Description: Calculate the sum-0f-squared haplotype frequency differences and mean per-site heterozygosity differences between the evolved populations and the base population.  

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: preprocess_step4-Splitting_haplotype_files.R <haplotypes> <founders> <helper>", call.=FALSE)
}

haplotypes <-  args[1]
founders <- args[2]
helper <- args[3]

# ============================================================================
# Load packages and sourced files

library(tictoc)
library(data.table)
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

print(projectDir)

indHapDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "/Individual_tx_tables/", analysisType = "hap_freqs", samplePattern = "^SEE12B02|^BAS02")

# ============================================================================
# Calculate the sum of squared differences between the haplotype frequencies of all samples and bas02 as well as average per-site heterozygosity differences

popDT <- do.call(rbind, indHapDTs)
popDT$population <- unlist(lapply(popDT$Chemical, function(x) strsplit(x, "_")[[1]][1]))
baseDTs <- popDT[id %like% "BAS02$"]
treatDTs <- popDT[!id %like% "BAS02$"]
wayDT <- split(treatDTs, treatDTs$population)
calcDifs <- lapply(wayDT, function(way) {
    base <- baseDTs[population == way$population[1]]
    weekDT <- split(way, way$Week)
    calcDifsLooper <- lapply(weekDT, function(week) {
        sampleSplit <- split(week, week$id)
        calculateHapDifsDTs <- Map(calcstats_ind.hap.dif.heterozygosities, sampleSplit, rep.list(base, sampleSplit), rep.list("id", sampleSplit))
        savingHapDifsDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "/Tables/Ind_tx_ind_haps_het_diffs_tables/", analysisType = "haps_het_diffs", idCol = "id")
        writeHapDifsDTs <- savingHapDifsDTs(calculateHapDifsDTs)
    } )
} )

# ============================================================================
# Trouble-shooting