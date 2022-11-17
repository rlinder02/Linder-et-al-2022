#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Finding aneuploidies
# Description: Find aneuploidies and infers their bounds using a Hidden Markov Model.

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
    stop("Usage: analysis_step14b-Aneuploidies_Summary_table.R <aneuploidies> <reps_using> <treatments> <treatment_key> <low_cov> <helper> <offsets>", call.=FALSE)
}

aneuploidies <-  args[1]
reps_using <- args[2]
treatments <- args[3]
treatment_key <- args[4]
low_cov <- args[5]
helper <- args[6]
offsets <- args[7]


# ============================================================================
# Load packages and sourced files

library(data.table)
library(tictoc)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()
options(digits = 10)

# ============================================================================
# Custom functions

integer_breaks <- function(n = 5, ...) {
    fxn <- function(x) {
        breaks <- floor(pretty(x, n, ...))
        names(breaks) <- attr(breaks, "labels")
        breaks
    }
    return(fxn)
}

# ============================================================================
# Load data

anDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "/", analysisType = "aneuploid_windows", samplePattern = "^SEE12B02.*[0-9]_aneuploid_windows_DT.txt$")
samplesUsing<- fread(reps_using, header = T)
treatmentDT <- fread(treatments, header = F)
treatKeyDT <- fread(treatment_key, header = T)
lowCovDTs <- fread(low_cov, header = T)
treatList <- treatmentDT[!V1 %in% lowCovDTs$id]
listOfTreatments <- treatList[V1 %like% "^SEE"][,V1]
treatKeyDT$Reps <- unlist(lapply(treatKeyDT$Abr, function(x) {
    length(grep(x, listOfTreatments))
} ) )
treatKeyDT$week <- substr(treatKeyDT$Abr, 1, 2)
treatKeyDT <- treatKeyDT[!Abr %like% "^[0-9][0-9]FA|^02NANA"]
offsets <- read.table(offsets, header = T)
offsets$lines <- ceiling(offsets[,2]/50)
offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
offsets$chr <- 1:17
g_l <- c(0, cumsum(offsets$len))
ch.bounds <- c(0, g_l[1:17] + offsets[,2])
ch.sizes <- diff(ch.bounds)
#correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
max.pos <- cumsum(ch.sizes)
mid.ch <- diff(ch.bounds)/2
midpt.ch <- ch.bounds[2:18] - mid.ch

# ============================================================================
#  Assemble a table of aneuploidies from tables generated for each replicate separately for each treatment
## classifying aneuploidies as whole or partial

allAns <- do.call(rbind, anDTs)
allAns$Pos <- trimws(format(allAns$POS, big.mark = ","))
allAns$avgNormCov <- unlist(lapply(allAns$avgNormCov, round, digits = 1))
repsMissing <- treatmentDT[which(!V1 %in% allAns$id)][V1 != "BAS02"][V1 != "DIP02"][!V1 %like% "YP000"][,V1]
repsMissing <- sort(rep(repsMissing, 2))
chrsMissing <- which(!(1:17) %in% allAns$chr)
chrsMissing <- sort(rep(chrsMissing, 2))
allAnsDT <- rbind(allAns, allAns[(.N+1):(.N+length(repsMissing))])
allAnsDT[which(is.na(id))] <- 0
allAnsDT[which(id == 0), "id"] <- unlist(lapply(repsMissing, function(x) {return(x)}))
extraChrRows <- length(repsMissing) %% length(chrsMissing)
allAnsDT[which(chr == 0), "chr"] <- unlist(lapply(c(rep(chrsMissing,length(repsMissing)/length(chrsMissing)), chrsMissing[seq_len(extraChrRows)]), function(x) {return(x)}))

allAnsDT$chemical <- unlist(lapply(1:nrow(allAnsDT), function(x){
    chemName <- paste(strsplit(allAnsDT$id[x], "")[[1]][9:10], collapse = "")
    return(chemName)
} ) )
anClassifier <- unlist(lapply(seq(1, nrow(allAnsDT), 2), function(x) {
    print(x)
    flush.console()
    size <- allAnsDT$POS[x+1] - allAnsDT$POS[x]
    chrSize <- ch.sizes[allAnsDT$chr[x]]
    relSize <- size/chrSize
    if(relSize >= .8) {return(rep("whole", 2))} else if(relSize == 0) {return(rep(0, 2))} else{return(rep("partial", 2))}
} ) )
allAnsDT$Classifier <- anClassifier
posRange <- unlist(lapply(seq(1, nrow(allAnsDT), 2), function(x) {paste0(allAnsDT$Pos[x], "-", allAnsDT$Pos[x+1])}))
posRange[posRange == "0-0"] <- 0
allAnsDT <- allAnsDT[seq(1, nrow(allAnsDT), 2)]    
allAnsDT$Pos <- posRange
allAnsDT$Class <- paste0(allAnsDT$Classifier, " ", allAnsDT$Type)
allAnsDTcpy <- copy(allAnsDT)
allAnsDT <- allAnsDT[chr != 17][!chemical %like% "FA[0-9][0-9]P"]
ansDT <- allAnsDT[, Replicate := lapply(.SD, function(x) {paste(strsplit(x, "")[[1]][c(14:16)], collapse = "")}), .SDcols = "id", by = seq_len(nrow(allAnsDT))
                  ][, week := lapply(.SD, function(x) {paste(strsplit(x, "")[[1]][c(4:5)], collapse = "")}), .SDcols = "id", by = seq_len(nrow(allAnsDT))
                    ][, chemical := lapply(.SD, function(x) {paste(strsplit(x, "")[[1]][4:10], collapse = "")}), .SDcols = "id", by = seq_len(nrow(allAnsDT)) 
                      ][, Treatment := lapply(.SD, function(x) {treatKeyDT[Abr == x, Chemical]}), .SDcols = "chemical", by = seq_len(nrow(allAnsDT))
                     ][, repFreq := paste0(Replicate, "(", avgNormCov, ")")
                      ][, chemId := substr(id, 9, 10)
                       ][, c("chr", "POS", "Class", "Replicate", "id", "Treatment", "Pos", "week", "repFreq", "chemId", "avgNormCov") 
                         ][order(Treatment, Replicate, chr, POS)]
ansDTdups <- ansDT[grepl("Duplication", Class) & week == 12 & id %in% samplesUsing$id]
ansDTdups[Class == "whole Duplication", Pos := 0]

ansDTdupsCpy <- copy(ansDTdups)
ansDTdupsCpy[, chrPos := paste0(chr, "_", Pos)]
maxCov <- ansDTdupsCpy[, max(avgNormCov), by = chrPos]
onecopydup <- maxCov[V1 >= 1.5, .N]
twocopydup <- maxCov[V1 >= 2, .N]

ansDTdups <- ansDTdups[, -c("avgNormCov")]
totalReps <- samplesUsing[, .N, by = chemical][, chemical := gsub("_12", "", chemical)]
formatDups <- ansDTdups[, toString(repFreq), by = .(Treatment, chr, Pos, Class, chemId)][order(chemId, chr, Pos)][, chr := as.roman(chr)]
setnames(formatDups, old = "V1", new = "Replicate(fold change)")
newDT <- formatDups[totalReps, on = c(Treatment = "chemical")][order(chemId)][,!c("Treatment")]
setnames(newDT, old = c("N", "chemId", "Pos"), new = c("Total_reps", "Chemical", "Interval"))
setcolorder(newDT, c("Chemical", "chr", "Interval", "Total_reps", "Replicate(fold change)"))
newDT2 <- newDT[Class == "whole Duplication", Interval := NA]
newDT2[, Class := gsub(" Duplication", "", Class)]
setnames(newDT2, old = c("Total_reps", "Class"), new = c("Total reps", "Duplication type"))
fwrite(newDT2, "Aneuploidy_table_v4.txt", sep = "\t", col.names = T)

# ============================================================================
# Trouble-shooting