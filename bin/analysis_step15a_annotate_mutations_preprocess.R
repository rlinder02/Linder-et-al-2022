#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Filtering mutations
# Description: Filter the list of potential mutations and convert into a format recognized by SNPeff.

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("Usage: analysis_step15a_annotate_mutations_preprocess.R <old_muts> <new_muts> <offsets> <helper>", call.=FALSE)
}

old_muts <-  args[1]
new_muts <- args[2]
offsets <- args[3]
helper <- args[4]

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

oldMuts <- fread(old_muts, header = T)
rhcMuts <- oldMuts[sample %like% "^RHC"]
baseMuts <- oldMuts[sample %like% "^postsp"]
newMuts <- fread(new_muts, header = T)
newMuts <- rbind(newMuts, rhcMuts, baseMuts)
offsets <- fread(offsets, header = T)
g_l <- posoff_chr.bounds(offsets)

# ============================================================================
# Filter and reformat mutations table into a snpeff compatible VCF file

newMuts <- newMuts[,chr2 := as.integer(as.roman(substr(chr, 4, nchar(chr))))
                   ][chr2 == 1000, chr2 := 17
                     ][order(chr2, pos, sample)]
newMuts <- newMuts[, gp := .(pos + g_l[chr2])]
newMuts <- newMuts[, id := sample][, sample := NULL]
newMuts$chemical <- unlist(lapply(1:nrow(newMuts), function(x){
    chemName <- paste(strsplit(newMuts$id[x], "")[[1]][9:10], collapse = "")
    return(chemName)
} ) )
newMuts <- newMuts[,mutID := paste(chr, pos, alt, sep = "_")]
newMuts <- newMuts[, c("fraction", "totalCov") := .(Nalt/(Nalt+Nref), Nref+Nalt)]
firstFilter <- newMuts[id %like% "^A[0-9]" | id %like% "^B[0-9]" | id %like% "^AB[0-9]" | id %like% "^BAS02$" | id %like% "^RHC" | id %like% "^DIP02$" | id %like% "^Base" | id %like% "FA02P" | id %like% "FA11P" | id %like% "^postsp"] # Filter out all mutations that show up in the original founders (including on plate controls - FA02 and FA11 founders here) as well as in any of the base populations, including clones directly derived from the base populations, as these mutations are likely misbehaving SNPs or were missed when sequencing the original founders or occur as part of the normal culturing of cells and sample prep for sequencing; are most interested in mutations that occur for specific treatments
newMuts <- newMuts[!mutID %in% unique(firstFilter$mutID)]
newMuts <- newMuts[totalCov >= 10] # Need at least 10x coverage at the mutations' position
newMuts <- newMuts[fraction >= 0.2] # At least 20% of the calls have to be for the mutation
countMuts <- newMuts[, .N, by = c("mutID", "chemical")]
uniqueMuts <- countMuts[, uniqueN(chemical), by = mutID]
goodMuts <- uniqueMuts[V1 <= 5] # Mutations that show up in more than five different treatments are unlikely to be specific for the treatment- more likely an artifact of the culturing conditions
goodMutsDT <- newMuts[mutID %in% goodMuts$mutID]
treatMutsCp <- copy(goodMutsDT)
treatMutsCp[chr == 'chrM', chr := 'chrmt']
treatMutsCp[, AllID := paste(chr, pos, id, ref, alt, sep = "_")]
goodMutsDT[, c("QUAL", "FILTER", "CHROM", "POS", "ID", "REF", "ALT", "INFO") := .(rep('.', .N), rep('.', .N), chr, pos, id, ref, alt, rep('.', .N))]
goodMutsDT[, c("Nref", "Nalt", "chr", "pos", "id", "ref", "alt", "mutID", "fraction", "chemical", "totalCov", "chr2", "gp") := .(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)]
setcolorder(goodMutsDT, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
goodMutsDT[CHROM == "chrM", CHROM := "chrMito"]
setnames(goodMutsDT, old = "CHROM", new = "#CHROM")
goodMutsDT <- goodMutsDT[-.N] # gets rid of the last empty row
goodMutsDF <- as.data.frame(goodMutsDT)
write.table(goodMutsDF, "SNP_mutations_sexual_restructured.vcf", row.names = FALSE, col.names = TRUE, quote = F, sep = "\t")
# ============================================================================
# Trouble-shooting