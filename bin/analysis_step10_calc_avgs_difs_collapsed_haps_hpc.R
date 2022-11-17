#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Averaging haplotype frequencies
# Description: Average haplotype frequencies across replicates within a chemical treatment; the sd is also calculated.

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Usage: analysis_step10_calc_avgs_difs_collapsed_haps_hpc.R <collapsed_haplotypes> ", call.=FALSE)
}

haplotypes <- args[1]

# ============================================================================
# Load packages and sourced files

library(data.table)
library(tidyverse)

# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

group.reps <- function(grouping_factor) {
    function(x) {
        groups <- split(x, grouping_factor)
        groups
    }
}

averaging.groups.base.difs <- function(list_of_reps, baseDT) {
    founder_avgs <- lapply(list_of_reps, function(x) {
        baseFreqs <- baseDT[gp == x$gp[1], "collapsedFreqs"]
        names(baseFreqs) <- "baseFreqs"
        bfreqs <- baseFreqs[, tstrsplit(baseFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
        freqs <- x[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
        founders <- as.character(x[1, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
        avgs <- freqs[, lapply(.SD, mean), .SDcols = c(paste0("V", 1:ncol(freqs)))]
        standDevs <- freqs[, lapply(.SD, sd), .SDcols = c(paste0("V", 1:ncol(freqs)))]
        #print(x$gp[1])
        #flush.console()
        difCols <- calculating.avg.difs(avgs, bfreqs)
        avgsCol <- unite(avgs, col = "avgFreq", sep = ";", remove = TRUE, na.rm = TRUE)
        sdCol <- unite(standDevs, col = "sdFreq", sep = ";", remove = TRUE, na.rm = TRUE)
        baseFreqs <- unite(bfreqs, col = "baseFreqs", sep = ";", remove = TRUE, na.rm = TRUE)
        info <- as.data.table(c(x[1,c("chr", "pos", "gp", "Chemical", "Replicate", "id", "Week", "chemWeek")], avgsCol, sdCol, baseFreqs, x[1, "collapsedFounders"]))
        infoAll <- cbind(info, difCols)
        infoAll[, avgSD := mean(as.numeric(standDevs))]
        infoAll
    } )
    concat <- do.call(rbind,founder_avgs)
    concat                    
}

calculating.avg.difs <- function(avgFreqs, baseFreqs) {
    avgDifs <- avgFreqs - baseFreqs
    avgDifsSq <- avgDifs^2
    sumSqDifs <- sum(avgDifsSq)
    names(sumSqDifs) <- "sumSqDifs"
    difSign <- avgDifs[, lapply(.SD, function(x) gsub("^-.*", "Decreasing", x))
    ][, lapply(.SD, function(x) gsub("^0$", "Neutral", x))
    ][, lapply(.SD, function(x) gsub("^([1-9].*|0.0*\\d+|0.\\d+)", "Increasing", x))]
    avgDifsCol <- unite(avgDifs, col = "avgDifs", sep = ";", remove = TRUE, na.rm = TRUE)
    avgDifsSqCol <- unite(avgDifsSq, col = "avgDifsSq", sep = ";", remove = TRUE, na.rm = TRUE)
    difSignCol <- unite(difSign, col = "difSign", sep = ";", remove = TRUE, na.rm = TRUE)
    difsDT <- as.data.table(c(avgDifsCol, avgDifsSqCol, difSignCol, sumSqDifs))
    difsDT
}


avg.sd.difs.replicates <- function(repsDT) {
    tictoc::tic()
    print(repsDT$Chemical[1])
    flush.console()
    chemDT <- repsDT[id != "BAS02"]
    baseDT <- repsDT[id == "BAS02"]
    reps_found <- unique(chemDT$id)
    if(length(reps_found) > 1) {
        cat('Calculating average and sd \n')
        flush.console()
        group.gp <- group.reps(chemDT$gp)
        list_of_reps <- group.gp(chemDT)
        avgDifDT <- averaging.groups.base.difs(list_of_reps, baseDT)
    } else{
        cat('Only one replicate \n')
        flush.console()
        avgDifDT <- repsDT
    }
    tictoc::toc()
    avgDifDT
}

# ============================================================================
# Load data needed for downstream analyses

hapFreqs <- fread(haplotypes, header = T, sep = "\t")

#hapFreqs <- read.table(hapDT, header = T, sep = "\t", as.is = T)
#hapFreqs <- as.data.table(hapFreqs)

# ============================================================================
# Normalize the collapsed founders across replicates within treatments using cutree as the key 

hapFreqsFilt <- hapFreqs
avgSDdifsDT <- avg.sd.difs.replicates(hapFreqsFilt)
fwrite(avgSDdifsDT, file=paste0(projectDir, "/Tables/Averaged_and_sd_tx_tables/", avgSDdifsDT$chemWeek[1], "_hap_freqs_avg_sd_difs_DT.txt"), col.names = T, sep = "\t")