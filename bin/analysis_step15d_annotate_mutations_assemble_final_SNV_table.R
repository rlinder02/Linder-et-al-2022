#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Annotate mutations
# Description: Further annotate the list of potential mutations using the S. cerevisiae gff file

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    stop("Usage: analysis_step15d_annotate_mutations_assemble_final_SNV_table.R <mutations_ann> <helper> <treatments> <treatment_key> <low_cov> <reps_using> ", call.=FALSE)
}

mutations <-  args[1]
helper <- args[2]
treatments <- args[3]
treatment_key <- args[4]
low_cov <- args[5]
reps_using <- args[6]

# ============================================================================
# Load packages and sourced files

library(tictoc)
library(dplyr)
library(data.table)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

# ============================================================================
# Load data

mutsDTcomp2 <- fread(mutations, header = T)
treatmentDT <- fread(treatments, header = F)
treatKeyDT <- fread(treatment_key, header = T)
lowCovDTs <- fread(low_cov, header = T)
treatList <- treatmentDT[!V1 %in% lowCovDTs$id]
listOfTreatments <- treatList[V1 %like% "^SEE"][,V1]
treatKeyDT$Reps <- unlist(lapply(treatKeyDT$Abr, function(x) {
    length(grep(x, listOfTreatments))
} ) )
treatKeyDT$week <- substr(treatKeyDT$Abr, 1, 2)
samplesUsing<- fread(reps_using, header = T)

# ============================================================================
# Assemble a table of SNVs for each replicate separately for each treatment; use gt table, which makes an html-formatted table with tons of customization
#filterYPDmuts <- mutsDTcomp2[chemical == "YP", "mutID"]
#mutsDTcomp2 <- mutsDTcomp2[!mutID %in% filterYPDmuts$mutID]
mutsDTcomp2$Pos <- format(as.numeric(mutsDTcomp2$POS), big.mark = ",")
mutsDTcomp2$fraction <- unlist(lapply(mutsDTcomp2$fraction, round, digits = 1))
mutsDTcomp2 <- mutsDTcomp2[, chrom := as.integer(as.roman(substr(Chr, 4, nchar(Chr))))][chrom == 1000, chrom := 17]
mutsDTcomp2$TF <- unlist(lapply(mutsDTcomp2$INFO, function(x) {
    tf <- strsplit(x, "\\|")[[1]][5]
    type <- strsplit(x, "\\|")[[1]][3]
    print(tf)
    flush.console()
    if(length(tf) > 0 & !is.na(tf) & type == "TFBS_Disruption") {return(tf)} else{
        return(NA)}
} ) )
mutsDTcomp2$chemical <- unlist(lapply(1:nrow(mutsDTcomp2), function(x){
    chemName <- paste(strsplit(mutsDTcomp2$ID[x], "")[[1]][9:10], collapse = "")
    return(chemName)
} ) )
mutsDTcpy <- copy(mutsDTcomp2)

mutDT <- mutsDTcomp2[!gene_hit %like% "TEL"]
mutDT <- mutDT[tx %like% "^SEE"]
mutDT <- mutDT[ID %in% samplesUsing$id]
treatKeyDT <- treatKeyDT[!Abr %like% "^[0-9][0-9]FA|^02NANA"]
#fwrite(mutDT, "SEE01_SNVsAnnotatedOriginal.txt", quote = F, sep = "\t", col.names = T)

mutDT <- mutDT[, Replicate := lapply(.SD, function(x) {paste(strsplit(x, "")[[1]][c(14:16)], collapse = "")}), .SDcols = "ID", by = seq_len(nrow(mutDT))
                 ][, chemical := lapply(.SD, function(x) {paste(strsplit(x, "")[[1]][c(9:10)], collapse = "")
               }), .SDcols = "ID", by = seq_len(nrow(mutDT))
                  ][, chr := as.numeric(as.roman(gsub("chr", "", Chr)))
                   #][, pos := format(Pos, big.mark = ",") 
                    ][, mutation := paste0(REF, ">", ALT)
                    ][, repFreq := paste0(Replicate, "(", fraction, ")")
                     ][, week := lapply(.SD, function(x) {paste(strsplit(x, "")[[1]][c(4:5)], collapse = "")}), .SDcols = "ID", by = seq_len(nrow(mutDT))
                        ][, c("ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "AllID", "tx", "dose", "tx_gene", "dose_gene", "fraction", "totalCov", "gp", "mutID", "Chr", "GO", "PATH", "Replicate", "chrom") := NULL
                       ][order(chemical, chr, POS)
                        ][, id := paste0(chemical, chr, POS, mutation, repFreq)]

mutDT2 <- mutDT[!duplicated(id)][, c("id", "POS", "week") := NULL]
totalReps <- samplesUsing[, .N, by = chemical][, chemical := gsub("_12", "", chemical)]
setnames(totalReps, old = "chemical", new = "Chemical")
totalReps2 <- merge(totalReps, treatKeyDT, by = "Chemical", all.x = T)[week == 12][, Chemical := substr(Abr, 6, 7)][, c("Abr", "V1", "week", "Reps") := NULL]
formatMuts <- mutDT2[, toString(repFreq), by = .(chemical, chr, Pos, mut_type, mutation, gene_hit, TF)][order(chemical, chr, Pos)][, chr := as.roman(chr)]
setnames(formatMuts, old = c("V1", "chemical"), new = c("Replicate(frequency)", "Chemical"))
newDT <- formatMuts[totalReps2, on = "Chemical"][order(Chemical)]
setnames(newDT, old = c("N", "mut_type"), new = c("Total_reps", "mutation_type"))
setcolorder(newDT, c("Chemical", "chr", "Pos", "Total_reps", "Replicate(frequency)", "mutation", "mutation_type", "gene_hit", "TF"))

avgFreq <- unlist(lapply(mutDT2$repFreq, function(x) {strsplit(x, "\\(")[[1]][2] } ) )
avgFreq2 <- as.numeric(gsub("\\)", "", avgFreq))
avgFreq3 <- mean(avgFreq2)
highFreqs <- length(which(avgFreq2 >= 0.5))
higherFreqs <- length(which(avgFreq2 >= 0.75))
explore <- mutDT[, .N, by = gene_hit]
multHit <- explore[N > 1]
multHit2 <- mutDT2[gene_hit %in% multHit$gene_hit][order(gene_hit)]

fwrite(newDT, "SNVsAnnotated_TFs.txt", sep = "\t", col.names = T)
# ============================================================================
# Trouble-shooting