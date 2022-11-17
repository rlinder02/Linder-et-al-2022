#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Preparing tables needed to collapse founder haplotypes
# Description: Preparing tables that are needed to collapse founder haplotypes that cannot be distinguished across the genome.

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    stop("Usage: preprocess_step7a_normalizing_collapsed_founders.R <reps_using> <founders> <haplotypes> <treatments> <treatment_key> <helper>", call.=FALSE)
}

reps_using <-  args[1]
founders <- args[2]
haplotypes <- args[3]
treatments <- args[4]
treatment_key <- args[5]
helper <- args[6]

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
hetDT <- fread(reps_using, header = T)
hapDTs <- fread(haplotypes, header=TRUE)
treatKeyDT <- fread(treatment_key, header = T)
treatmentDT <- fread(treatments, header = F)
listOfTreatments <- treatmentDT[V1 %in% hetDT$id | V1 %like% "^BAS02$"][,V1]

# ============================================================================
# Prep files to collapse founders
            
treatIds <- data.table(gsub("R[0-9][0-9]$", "", listOfTreatments))[,outcome:= rleid(V1)]
listOfTreatmentsDT <- data.table(listOfTreatments)[, outcome := treatIds$outcome]
grpIdx <- split(listOfTreatmentsDT, listOfTreatmentsDT$outcome)
grpIdxBase <- lapply(grpIdx, function(x) {
  newList <- rbindlist(list(x, list("BAS02", "BAS02")))
  newList$listOfTreatments
} )
grpIdxBase <- grpIdxBase[-1]
chemsUsing <- unlist(lapply(grpIdxBase, function(x) {
  chem <- gsub("R[0-9][0-9]$", "", x[1])
  chem
} ) )
chemsDT <- as.data.table(chemsUsing)
newHapDT <- hapDTs[poolroot %in% listOfTreatmentsDT$listOfTreatments]
fwrite(chemsDT, "SEE01_chems_revised.txt", col.names = F)
fwrite(newHapDT, "SEE01_reps_using_haps.txt", col.names = T)

### Run the collapse_founder_haps.sub file on the hpc (first run "dos2unix" on file) which runs the SEE01_collapsing_haplotypes_hpc.R script; the output will be used as input for subsequent steps and should be transferred to the "All_reps_tx_tables" folder

# ============================================================================
# Trouble-shooting

