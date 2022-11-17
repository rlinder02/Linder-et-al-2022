#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Splittig haplotype files
# Description: Split the haplotype frequency calls into individual files for each sample

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("Usage: preprocess_step4-Splitting_haplotype_files.R <haplotypes> <low_coverage_samples> <founders> <helper>", call.=FALSE)
}

haplotypes <-  args[1]
low_cov <- args[2]
founders <- args[3]
helper <- args[4]

# ============================================================================
# Load packages and sourced files

library(tictoc)
library(data.table)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()
options(digits = 10)

# ============================================================================
# Custom functions

# ============================================================================
# Load data

# ============================================================================
# Run analyses

hapDTs <- fread(haplotypes, header=TRUE)
lowCovDT <- fread(low_cov, header = T)
lowCovDT[, id := gsub("_", "", id)]
founderNames <- fread(founders, header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)

# ============================================================================
# Filter out low coverage samples

hapDT <- hapDTs[!lowCovDT, on = .(poolroot = id)]
findNAs <- hapDT[founderfreqs %like% "NA", unique(poolroot)]
hapDT <- hapDT[!poolroot %in% findNAs]
samples <- as.data.table(unique(hapDT$poolroot))
samplesCpy <- copy(samples)
treatments <- samplesCpy[, c('V1') := lapply(V1, function(x) strsplit(x, "R[0-9][0-9]$")[[1]][1])][, .(V1 = unique(V1))]
fwrite(samples, "SEE01_unique_samples.txt", sep = '\t', col.names = F)
fwrite(treatments, "SEE01_unique_treatments.txt", sep = '\t', col.names = F)

# ============================================================================
# Make a look-up table of all chemical treatments and their abbreviations

treatsDT <- as.data.table(treatments)
print(treatsDT)
treatKeyDT <- treatsDT[, Abr := lapply(.SD, function(x) {
    x <- x[[1]]
    print(x)
    if(nchar(x) < 10) {x} else {
     paste(strsplit(x, "")[[1]][9:10], collapse = "")} } ), .SDcols = "V1", by = seq_len(nrow(treatsDT))
    ][, Chemical := unlist(lapply(seq_len(nrow(treatsDT)), function(x) {
        if(Abr[x] == "TU") {return("tunicamycin")} else if(Abr[x] == "CA") {return("caffeine")} else if(Abr[x] == "CD") {
            return("cadmium_chloride")} else if(Abr[x] == "CP") {return("chlorpromazine")} else if(Abr[x] == "CI") {
                return("cisplatin")} else if(Abr[x] == "FL") {return("fluconazole")} else if(Abr[x] == "DM") {
                    return("dmso")} else if(Abr[x] == "ET") {return("ethanol")} else if(Abr[x] == "GA") {
                        return("glacial_acetic_acid")} else if(Abr[x] == "NT") {
                            return("nicotine")} else if (Abr[x] == "NC") {
                                return("sodium_chloride")} else if(Abr[x] == "NM") {
                                    return("nicotinamide")} else if(Abr[x] == "SO") {
                                        return("sodium_sulfite")} else if(Abr[x] == "UR") {
                                            return("urea")} else if(Abr[x] == "DI") {
                                                return("diamide")} else if(Abr[x] == "YP") {
                                                    return("YPD")} else {
                                                        return(NA)} } ) ) ] 
treatKeyDT[V1 == "BAS02", c("Abr", "Chemical") := .("BAS02", "18way_base")][
    V1 == "DIP02", c("Abr", "Chemical") := .("DIP02", "2way_base")]
treatKeyDT$Chemical <- unlist(lapply(1:nrow(treatKeyDT), function(x) {
    if(grepl("B02", treatKeyDT$V1[x])) {
        return(paste0("18way_", treatKeyDT$Chemical[[x]]))} else if(grepl("D02", treatKeyDT$V1[x])) {
            return(paste0("2way_", treatKeyDT$Chemical[[x]]))} else {
                return(treatKeyDT$Chemical[[x]])} 
} ) )
print(treatKeyDT)
treatKeyDT$Abr <- unlist(lapply(treatKeyDT$V1, function(x) {
    paste(strsplit(x, "")[[1]][4:10], collapse = "")
}))

fwrite(treatKeyDT, "treatment_key.txt", col.names = T)
# ============================================================================
# Split haplotype calls into separate files for each individual sample

treatmentDT <- samples
listOfTreatments <- treatmentDT[V1 %like% "^SEE12B02" | V1 %like% "^BAS02$" | V1 %like% "^DIP02"][,V1]
splittingIndSamples <- filesplitter_split.collapsed.haps.basic.samples(listOfTreatments, idCol = "poolroot")
splitIndSamplesDTs <- splittingIndSamples(hapDT)
savingIndDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "/Tables/Individual_tx_tables/", analysisType = "hap_freqs", idCol = "id")
writeIndDTs <- savingIndDTs(splitIndSamplesDTs)

# ============================================================================
# Trouble-shooting