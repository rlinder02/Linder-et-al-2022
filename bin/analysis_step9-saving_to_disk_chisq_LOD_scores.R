#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Chi-Square test
# Description: Running a genome-wide chi-square test to find regions in evolved populations with haplotype frequencies significantly different from the base population that show evidence of selection

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: analysis_step9-saving_to_disk_chisq_LOD_scores.R <transformed_haplotypes> <founders> <helper>", call.=FALSE)
}

haplotypes <- args[1]
founders <- args[2]
helper <- args[3]

# ============================================================================
# Load packages and sourced files

library(tictoc)
library(data.table)
library(tidyverse)
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
indHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "/", analysisType = "haps_sq_diffs", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]_")

# ============================================================================
# Calculate the chi square LOD scores of samples vs bas02 

popDT <- do.call(rbind, indHapDifsDTs)
print(popDT)
flush.console()

wayDT <- split(popDT, popDT$population)
calcDifs <- lapply(wayDT, function(way) {
    weekDT <- split(way, way$Week)
    calcDifsLooper <- lapply(weekDT, function(week) {
        sampleSplit <- split(week, week$chemWeek)
        calculateHapDifsDTs <- Map(calcstats_asin.sqrt.hap.avg.dif.transforms, sampleSplit, rep.list("id", sampleSplit), rep.list(0.03, sampleSplit), rep.list(0.005, sampleSplit))
        savingHapDifsDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "/Tables/Hap_adjusted_LOD_tables/", analysisType = "hap_adjusted_LOD", idCol = "id")
        writeHapDifsDTs <- savingHapDifsDTs(calculateHapDifsDTs)
    } )
} )

# ============================================================================
# Calculate the chi square LOD scores of all chemicals split into two groups of three (YPD is two groups of two) to get at repeatability

wayDT <- split(popDT, popDT$population)
calcDifs <- lapply(wayDT, function(way) {
    weekDT <- split(way, way$Week)
    calcDifsLooper <- lapply(weekDT, function(week) {
        sampleSplit <- split(week, week$chemWeek)
        sampleSubSplit <- lapply(sampleSplit, function(samp) {
            reps <- samp[, unique(Replicate)]
            if(samp$Chemical[1] == "18way_YPD") {
                set1 <- reps[1:2]
                set2 <- reps[3:4]
            } else {
                set1 <- reps[1:3]
                set2 <- reps[4:6]
            }
            setLst <- list(set1, set2)            
            repLst <- lapply(1:length(setLst), function(x) {
                #print(x)
                #flush.console()
                newDT <- samp[Replicate %in% setLst[[x]]]
                newDT$Chemical <- paste0(newDT$Chemical[1], "_grp_", x)
                newDT
            } )
            calculateLODDTs <- Map(calcstats_asin.sqrt.hap.avg.dif.transforms, repLst, rep.list("id", repLst), rep.list(0.03, repLst), rep.list(0.005, repLst))
            savingLODDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "/Tables/Repeatability_tables/", analysisType = "hap_adjusted_LOD", idCol = "Chemical")
            writeLODDTs <- savingLODDTs(calculateLODDTs)
        } )
    } )
} )


# ============================================================================
# Calculate the chi square LOD scores of cadmium chloride split into two groups of five vs bas02 to show get similar results

cadDT <- popDT[Chemical == "18way_cadmium_chloride"]
if(cadDT[, .N] > 2) {
    cadReps <- split(cadDT, cadDT$Replicate)
    cadReps1 <- do.call(rbind, cadReps[1:5])
    cadReps2 <- do.call(rbind, cadReps[6:10])
    cadGrps <- list(cadReps1, cadReps2)
    counter <- 0
    calcDifsLooper <- lapply(cadGrps, function(cad) {
        counter <<- counter + 1
        calculateHapDifsDT <- calcstats_asin.sqrt.hap.avg.dif.transforms(cad, "id", 0.03, 0.005)
        fileName <- paste0("cad_pt", counter, ".txt")
        fwrite(calculateHapDifsDT, file = paste0(projectDir, "/Tables/Pleiotropy_tables/", fileName), sep = "\t", col.names = T)
    } )
} else {
    emptyData <- data.table(empty = 0)
    fileName <- "empty.txt"
    fwrite(emptyData, file = paste0(projectDir, "/Tables/Pleiotropy_tables/", fileName), sep = "\t", col.names = T)
}

# ============================================================================
# Trouble-shooting