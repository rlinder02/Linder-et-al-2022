#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Filtering replicate populations  
# Description: Making a look-up table of samples included in downstream analyses.  

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    stop("Usage: preprocess_step6-heterozygosity_filtering_and_further_analyses.R <hap_diffs> <founders> <treatments> <treatment_key> <helper>", call.=FALSE)
}

hap_diffs <-  args[1]
founders <- args[2]
treatments <- args[3]
treatment_key <- args[4]
helper <- args[5]

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

indHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "/Ind_tx_ind_haps_het_diffs_tables/", analysisType = "haps_het_diffs", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]_")
founderNames <- fread(founders, header = F)
founderNames <- founderNames[,V1]
treatmentDT <- fread(treatments, header = F)
treatKeyDT <- fread(treatment_key, header = T)

# ============================================================================
# Create a master list of chemicals and replicates using that includes heterozygosity information

popDT <- do.call(rbind, indHapDifsDTs)
popDT$treatment <- sub("_", "-", popDT$Chemical)
popDT$treatment <- gsub("^(.*)-", "", popDT$treatment)

baseDT <- popDT[chemWeek == "18way_caffeine_12" & Replicate == "R01"][, c("heterozygosity", "Replicate", "Chemical", "treatment") := .(baseHeterozygosityCol, "BAS02", "18way_base_2", "BAS02")]
baseDTCpy <- copy(baseDT)
baseDT2 <- baseDT[, c("chemical", "Replicate", "hetChange", "het", "id", "population", "Week") := .(Chemical, Replicate, 0, heterozygosity, treatment, strsplit(Chemical[1], "_")[[1]][1], 0)][, c("chemical", "Replicate", "hetChange", "het", "id", "population", "Week")][, hetUp := NA]

waySplit <- split(popDT, popDT$population)

wayLoop <- lapply(waySplit, function(wayDT) {
  #base <- baseDTs[population == way$population[1]]
  weekSplit <- split(wayDT, wayDT$Week)
  weekLoop <- lapply(weekSplit, function(weekDT) {
    chemSplit <- split(weekDT, weekDT$treatment)
    chemLooper <- lapply(chemSplit, function(chem) {
      print(chem$chemWeek[1])
      flush.console()
      avgHet <- chem[, mean(heterozygosity), by = id]
      baseHet <- baseDT2[, .(id = "BAS02", V1 = mean(het))]
      allHet <- rbind(baseHet, avgHet)
      allChem <- rbind(baseDTCpy, chem)
      uniqueChem <- unique(allChem$chemWeek)
      chemName <- ifelse(length(uniqueChem) > 1,  uniqueChem[2], uniqueChem)
      hetChange <- allChem[, unique(hetAvgDevCol), by = Replicate][2:.N, V1]
      hetUp <- allChem[, unique(hetUpCol), by = Replicate][2:.N, V1]
      if(chemName == "18way_caffeine_12") {
        id <- c("BAS02", unique(allChem$id))} else {
          id <- c("BAS02", unique(allChem$id)[-1])}
      hetChangeDT <- data.table(chemical = chemName, Replicate = unique(allChem$Replicate), hetChange = c(0, hetChange), het = allHet$V1, hetUp = c(NA, hetUp), id = id, population = allChem$population[1], Week = allChem$Week[1])
      hetChangeDT
    })
    weekBind <- do.call(rbind, chemLooper)
  })
  wayBind <- do.call(rbind, weekLoop)
  wayBind
})
allBind <- do.call(rbind, wayLoop)
allBind[1, chemical := "BAS02"]
base <- allBind[1]
allBind2 <- allBind[Replicate != "BAS02"]
allBind3 <- rbind(base,allBind2)
allBind3[, reason := "heterozygosity"]
allBind3[id %like% "SEE12B02CD600R01|SEE12B02CD600R02|SEE12B02CD600R03|SEE12B02CP010R12|SEE12B02CP010R13|SEE12B02GA120R08|SEE12B02YP000R10", reason := "correlation"]
topTen <- allBind3[chemical %like% "cadmium|chlorpromazine|diamid|sodium_chloride|glacial"][order(chemical, hetChange)][, head(.SD, 10), by = chemical]
topSix <- allBind3[chemical %like% "urea|YPD"][order(chemical, hetChange)][, head(.SD, 6), by = chemical]
repsUsing <- rbind(topTen, topSix)
repsUsing <- repsUsing[!id %like% "SEE12B02CD600R01|SEE12B02CD600R02|SEE12B02CD600R03|SEE12B02CP010R12|SEE12B02CP010R13|SEE12B02GA120R08|SEE12B02YP000R10"]
repsUsingFinal <- repsUsing[order(chemical, Replicate)][, -c("reason")]
fwrite(repsUsingFinal, paste0(projectDir, "/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"), quote = F, sep = "\t", col.names = T)
fwrite(allBind3, paste0(projectDir, "/Tables/Summary_tables/SEE01_All_Reps_Hets.txt"), quote = F, sep = "\t", col.names = T)

# ============================================================================
# Trouble-shooting