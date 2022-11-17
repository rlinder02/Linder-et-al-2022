#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2020
# Title:       Genomewide coverage
# Description: Calculates the mean genomewide coverage per sample, excluding the mitochondria and repetitive regions, and stores this as a table. Also makes separate plots for each sample. Plots are iterated through chromosome by chromosome in a specified window.

##############################################################################

# ============================================================================
# Load packages and sourced files

library(data.table)
library(R.utils)

# ============================================================================
# Set global options

options(digits = 10)


# ============================================================================
# Custom functions


# ============================================================================
# Load data

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: preprocess_step1-reformat_snp_files.r <snp_counts> <snp_frequencies> <offsets>", call.=FALSE)
}

counts <-  args[1]
freqs <- args[2]
offsets <- args[3]
snp_freqs <- fread(freqs, header=TRUE)
snp_counts <- fread(counts, header=TRUE)
offsets <- read.table(offsets, header = T)
offsets$lines <- ceiling(offsets[,2]/50)
offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
offsets$chr <- 1:17
g_l <- c(0, cumsum(offsets$len))
ch.bounds <- c(0, g_l[1:17] + offsets[,2])
#correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
max.pos <- g_l[length(g_l)]

# ============================================================================
# Run analyses
## First calculate coverages genomewide and write the tables to disk.


#snp_data <- snp_data[,-c(grep("_AEE", names(snp_data))), with = FALSE]
snp_counts <- snp_counts[,chr := as.integer(as.roman(substr(CHROM, 4, nchar(CHROM))))
    ][chr == 1000, chr := 17
     ][,CHROM := NULL
      ][order(chr, POS)
        ][, gp := POS + g_l[chr]
         ]
snp_freqs <- snp_freqs[,chr := as.integer(as.roman(substr(CHROM, 4, nchar(CHROM))))
    ][chr == 1000, chr := 17
     ][,CHROM := NULL
      ][order(chr, POS)
       ][, gp := POS + g_l[chr]
        ]

#call_nums <- grep("^N_", names(snp_data))
#snp_freqs <- grep("^freq_", names(snp_data))
#snp_coverage <- snp_data[,-c(2,3,snp_freqs), with = FALSE]
re_ordered_cols <- c("chr", "POS", "gp")
#snp_data <- snp_data[,-c(call_nums), with = FALSE
#		][,REF := as.character(REF)
#		 ][,ALT := as.character(ALT)
#		  ]
#setcolorder(snp_data, c(re_ordered_cols, setdiff(names(snp_data), re_ordered_cols)))
#setcolorder(snp_coverage, c(re_ordered_cols, setdiff(names(snp_coverage), re_ordered_cols)))
setcolorder(snp_counts, c(re_ordered_cols, setdiff(names(snp_counts), re_ordered_cols)))
setcolorder(snp_freqs, c(re_ordered_cols, setdiff(names(snp_freqs), re_ordered_cols)))

ind_snp_data <- copy(snp_freqs)
#name_pruning <- grep("^freq", names(ind_snp_data))
#cov_name_pruning <- grep("^N", names(snp_coverage))
#new_names <- as.character(vapply(names(ind_snp_data)[name_pruning], function(x) strsplit(x, "_")[[1]][2], character(1)))
#new_cov_names <- as.character(vapply(names(snp_coverage)[cov_name_pruning], function(x) strsplit(x, "_")[[1]][2], character(1)))
#names(ind_snp_data)[name_pruning] <- new_names
#names(snp_coverage)[cov_name_pruning] <- new_cov_names
snp_freqs <- snp_freqs[, !grepl("^o", names(snp_freqs)), with = FALSE]
find_replicates <- snp_freqs[,grepl("R0[1-9]$", names(snp_freqs)), with = FALSE]
replicate_indices <- grep("R0[1-9]$|R1[0-9]$", names(snp_freqs))
names(find_replicates) <- substr(names(find_replicates), 1, (nchar(names(find_replicates)) - 3))
groups <- match(names(find_replicates), names(find_replicates))
find_replicates <- rbind(find_replicates, as.list(c(groups)))
find_replicates <- as.data.frame(find_replicates, stringsAsFactors = FALSE)
snp_means <- lapply(unique(groups), function(x) {
    extract_columns <- which(find_replicates[nrow(find_replicates),] == x)
    avgs <- rowMeans(find_replicates[extract_columns], na.rm = TRUE)
    names(avgs) <- rep(names(find_replicates)[x], length(avgs))
    avgs 
} ) 
values <- Map(as.numeric, snp_means)
avg_df <- as.data.frame(values)
names(avg_df) <- unlist(lapply(snp_means, function(x) {names(x)[1]}))
avg_df <- avg_df[-nrow(avg_df),]
snp_freqs <- as.data.frame(snp_freqs, stringsAsFactors = FALSE)
snp_freqs <- snp_freqs[,-c(replicate_indices)]
#name_pruning2 <- grep("^freq", names(snp_data))
#new_names2 <- as.character(vapply(names(snp_freqs)[name_pruning2], function(x) strsplit(x, "_")[[1]][2], character(1)))
#names(snp_data)[name_pruning2] <- new_names2
avg_snp_data <- cbind(snp_freqs, avg_df)

fwrite(snp_counts,"SNPtable.Dec28.sort_coverage_restructured.txt", quote = F, sep = "\t", col.names = T, row.names = F)

fwrite(ind_snp_data,"SNPtable.Dec28.sort_ind_restructured.txt", quote = F, sep = "\t", col.names = T, row.names = F)

fwrite(avg_snp_data,"SNPtable.Dec28.sort_avg_restructured.txt", quote = F, sep = "\t", col.names = T, row.names = F)