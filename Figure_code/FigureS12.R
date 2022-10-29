##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2020
# Title:       Genomewide coverage
# Description: Calculates the mean genomewide coverage per sample, excluding the mitochondria and repetitive regions, and stores this as a table. Also makes separate plots for each sample. Plots are iterated through chromosome by chromosome in a specified window.

##############################################################################

# ============================================================================
# Load packages and sourced files
## Sourced files are kept in the default working directory	

library(ggplot2)
library(data.table)
library(tictoc)
library(cowplot)
source('plotting/Plotting_genomewide_statistics.R')
source('seq/Position_offsetter.R')
source('seq/Coverage_calculator.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')

# ============================================================================
# Set global options
defDir <- getwd()
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

options(digits = 10)

# ============================================================================
# Custom functions

pattern_specifier <- function(pattern) {
    function(x) {
        x[grepl(pattern, x)]
    }
}

reading_in_dfs_as_a_list <- function(directory, file_prefix, pattern, file_suffix) {
    if(substr(getwd(), 1, 3) == "C:/") {
        setwd(file.path("C:/Users/rob/Dropbox/Long_lab/Papers/Resource_paper/Results_data/", directory)) } else {
            setwd(file.path("/Users/robertlinder/Dropbox/Long_lab/Papers/Resource_paper/Results_data/", directory)) 
        }
    sample_pattern <- pattern_specifier(pattern)
    samples <- sample_pattern(dir())
    read_in_files <- lapply(samples, function(x) {
        print(x)
        flush.console()
        data <- read.table(x, header = T, sep = "\t", check.names = FALSE)
        i <- sapply(data, is.factor)
        data[i] <- lapply(data[i], as.character)
        data
    } )
    if(substr(getwd(), 1, 3) == "C:/") {
        setwd("C:/Users/rob/Dropbox/Long_lab/Papers/Resource_paper/Results_data/Sequencing_data/Processed_sequencing_data") } else {
            setwd("/Users/robertlinder/Dropbox/Long_lab/Papers/Resource_paper/Results_data/Sequencing_data/Processed_sequencing_data")
        }
    read_in_files
}

# ============================================================================
# Run analyses
setwd("/Users/robertlinder/Dropbox/Long_lab/Papers/Resource_paper/Results_data/Sequencing_data/Processed_sequencing_data")
bas02PrivSnps <- reading_in_dfs_as_a_list(directory = 'Analyzed_results_data', pattern = 'BAS02.txt$')
bas02PrivSnpsDF <- bas02PrivSnps[[1]]
bas02PrivSnpsDT <- as.data.table(bas02PrivSnps)
bas02PrivSnpsDT2 <- bas02PrivSnpsDT[unique_snp_freq > 0]
bas02LowFreqsDT <- bas02PrivSnpsDT[unique_snp_freq <= 0.04] 
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

avgMAF <- bas02PrivSnpsDT2[, mean(unique_snp_freq)]
maf004 <- bas02PrivSnpsDT2[unique_snp_freq > 0.004, .N]/bas02PrivSnpsDT2[, .N]
maf0004 <- bas02PrivSnpsDT2[unique_snp_freq > 0.0004, .N]/bas02PrivSnpsDT2[, .N]

bas02PrivSnpsDT[, unique_snp_count := unique_snp_freq*bas02_snp_cov]
# ============================================================================
# Histogram of coverages genomewide

mafPlotA <- ggplot() + geom_histogram(data=bas02PrivSnpsDT[unique_snp_freq <= 0.5], aes(x = unique_snp_freq,  y = ..count../(sum(..count..))), binwidth = 0.049, center = 0.025) + scale_x_continuous(breaks=seq(0.001, 0.5, by = 0.049)) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + xlab("MAF") + ylab("Density")

mafPlotB <- ggplot() + geom_histogram(data=bas02LowFreqsDT, aes(x = unique_snp_freq, y = ..count../(sum(..count..))), binwidth = 0.003, center = 0.0025) + scale_x_continuous(breaks=seq(0.001, max(bas02LowFreqsDT$unique_snp_freq) + 0.003, by = 0.003)) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + xlab("MAF") + ylab("Density")


savePlot <- plot_grid(mafPlotA, mafPlotB, align = 'vh', labels = c("A", "B"), ncol = 1, hjust = -1)

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureSx_MAF_private_snps_hist.pdf"), savePlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureSx_MAF_private_snps_hist.png"), savePlot, width = 8, height = 9, units = "in", dpi = 350) 

visibleSites <- bas02PrivSnpsDT[unique_snp_freq >= 0.004, .N]/bas02PrivSnpsDT[, .N]

# ============================================================================
# Histogram of minor allele counts for private SNPs vs counts (so for minor allele counts from 0-whatever, how many times are there 0, 1, 2, etc... counts of each minor allele?)

mafPlotA <- ggplot() + geom_histogram(data=bas02PrivSnpsDT[chr != 17], aes(x = unique_snp_count,  y = ..count..), binwidth = 1) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + xlab("Minor allele count") + ylab("Count") + scale_x_continuous(limits = c(-1, 11), breaks = seq(0,10,2), labels = seq(0,10,2))

mafPlotB <- ggplot() + geom_histogram(data=bas02PrivSnpsDT[chr != 17 & unique_snp_count >= 10], aes(x = unique_snp_count,  y = ..count..), binwidth = 1) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + xlab("Minor allele count") + ylab("Count") + scale_x_continuous(limits = c(10, 1000), breaks = c(10, 200, 400, 600, 800, 1000), labels = c(10, 200, 400, 600, 800, 1000))


savePlot <- plot_grid(mafPlotA, mafPlotB, align = 'vh', labels = c("A", "B"), ncol = 1, hjust = -1)

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureSx_MAF_private_snps_hist_counts.pdf"), savePlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureSx_MAF_private_snps_hist_counts.png"), savePlot, width = 8, height = 9, units = "in", dpi = 350) 

visibleSites <- bas02PrivSnpsDT[unique_snp_freq >= 0.004, .N]/bas02PrivSnpsDT[, .N]
noCov <- bas02PrivSnpsDT[chr != 17 & unique_snp_freq == 0, .N]/bas02PrivSnpsDT[chr != 17, .N]
# ============================================================================
# Trouble-shooting
