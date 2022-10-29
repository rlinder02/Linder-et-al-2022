#####################  Loading required packages  ############################################

library(data.table)
library(viridis)
library(tictoc)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(readxl)
library(gprofiler2)
library(gridExtra)
library(egg)
library(scales)
library(dplyr)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggcorrplot)
library(ggbeeswarm)
library(tidyverse)
library(DescTools)
library(lemon)
library(abind)
library(cowplot)
library(tictoc)
library(hexbin)
library(ggpmisc)

source('utils/Reading_files_helper.R')

#######################  Setting global options  ################################################

options(digits = 10)
mypal <- rep(brewer.pal(8,'Dark2'),3)
defDir <- getwd()
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

#############  Custom Functions  #################################################################

remove_factors <- function(df) {
	x <- sapply(df, is.factor)
	df[x] <- lapply(df[x], as.character)
	df
}

renaming_founders <- function(old_names, new_names) {
	function(x) {
		names(x)[names(x) %in% old_names] <- new_names
		x
	}	
}

repeating_list <- function(object, times_to_repeat) {
	rep(list(object), times_to_repeat)
}

pattern_specifier <- function(pattern) {
	function(x) {
		x[grepl(pattern, x)]
	}
}

private_snp_freqs_vs_hap_freqs <- function(treatment, sample, snp_df, snp_cov_df, private_snp_df, merged_founders, collapsed_founders, bottom_cov_cutoff, top_cov_cutoff) {
    tic('time per sample')
    hapDF <- sample
    print(hapDF$treatment[1])
    flush.console()
    snp_metadata <- snp_df[,1:3]
    snp_cov <- snp_cov_df[,grep(treatment, names(snp_cov_df))]
    snp_data <- snp_df[,grep(treatment, names(snp_df))]
    snp_df2 <- cbind(snp_metadata, snp_data, snp_cov)
    private_snps <- snp_df2[snp_df2$gp %in% private_snp_df$gp,]
    private_snps_overlap <- private_snp_df[private_snp_df$gp %in% private_snps$gp,]
    private_snps <- cbind(private_snps, private_snps_overlap[,c('unique_snp', 'unique_founder')])
    private_snps <- private_snps[!private_snps$unique_founder %in% merged_founders,]
    #bottom_cov_cutoff1 <- quantile(bas02_private_snps$bas02_snp_cov, bottom_cov_cutoff, na.rm = T)
    #top_cov_cutoff <- quantile(bas02_private_snps$bas02_snp_cov, top_cov_cutoff)
    private_snps <- private_snps[private_snps$snp_cov >= bottom_cov_cutoff,] 
    #& bas02_private_snps$bas02_snp_cov <= top_cov_cutoff,]
    private_snps$unique_snp_freq <- 0
    private_snps$unique_snp_freq[private_snps$unique_snp == 0] <- 1 - private_snps$snp_data[private_snps$unique_snp == 0]
    private_snps$unique_snp_freq[private_snps$unique_snp == 1] <- private_snps$snp_data[private_snps$unique_snp == 1]
    private_snps$hap_calls <- unlist(lapply(1:nrow(private_snps), function(x) {
        #print(x)
        left_of <- which(hapDF$gp < private_snps$gp[x])
        right_of <- which(hapDF$gp > private_snps$gp[x])
        if(length(left_of) > 0 && length(right_of) > 0) {
            less_than <- max(left_of)
            more_than <- min(right_of)
            if(private_snps$chr[x] == hapDF$chr[less_than]) {
                hap_row <- hapDF[less_than,]} else {
                    hap_row <- hapDF[more_than,] } 
        } else if(length(left_of) < 1 && length(right_of) > 0) {
            more_than <- min(right_of)
            hap_row <- hapDF[more_than,] } else if(length(left_of) > 0 && length(right_of) < 1) {
                less_than <- max(left_of)
                hap_row <- hapDF[less_than,] } 
        hap_idx <- which(names(hap_row) %in% private_snps$unique_founder[x])
        hap_row <- as.data.frame(hap_row)
        hap_value <- hap_row[1,hap_idx]
        hap_value
    } ) )
    private_snps$bins <- 1
    private_snps$bins[private_snps$unique_founder %in% collapsed_founders] <- 2
    private_snps$difs <- private_snps$hap_calls - private_snps$unique_snp_freq
    toc()
    private_snps
}

find_snps_to_filter_out <- function(window, error_rate, private_snp_df, bas02_private_snp_df) {
	errs_df <- bas02_private_snp_df[bas02_private_snp_df$bins == 1 & bas02_private_snp_df$difs > error_rate | bas02_private_snp_df$difs < -(error_rate),]
	bad_snps <- unlist(lapply(1:nrow(errs_df), function(x) {
		#print(x)
		chrom <- errs_df$chr[x]
		gposition <- errs_df$gp[x]
		founder <- errs_df$unique_founder[x]
		hap_freq <- errs_df$hap_calls[x]
		snp_freq <- errs_df$unique_snp_freq[x]
		lower_bound <- gposition - window
		upper_bound <- gposition + window
		raw_interval <- private_snp_df[private_snp_df$gp >= lower_bound & private_snp_df$gp <= upper_bound,]
		interval <- raw_interval[raw_interval$chr == chrom,]
		nearby_snps <- interval[interval$unique_founder == founder & interval$gp != gposition,]
		if(nrow(nearby_snps) > 0) {
			nearby_snps_haps <- bas02_private_snp_df[bas02_private_snp_df$gp %in% nearby_snps$gp,]
			if(nrow(nearby_snps_haps) > 0) {
				nearby_snps_freqs <- mean(nearby_snps_haps$unique_snp_freq)
				if(abs(hap_freq - nearby_snps_freqs) < 0.01 || abs(snp_freq - nearby_snps_freqs) > 0.05) {return(gposition)}
			}	
		}
	} ) )
	bad_snps
}

delete_snps <- function(snp_list) {
	function(x) {
		if(length(snp_list) > 0) {
		df_corrected <- x[!x$gp %in% snp_list,] } else {
			df_corrected <- x}
	df_corrected
	}
}

saving_lists_of_dfs_to_txt_files <- function(list_of_dfs, file_prefix, sample_names, file_suffix) {
    counter <- 0
    lapply(list_of_dfs, function(x) {
        counter <<- counter + 1
        print(sample_names[counter])
        names_of_files <- paste0(file_prefix, sample_names[counter], ".txt")
        print(names_of_files)
        flush.console()
        write.table(x, names_of_files, sep = '\t', col.names = T, row.names = F, quote = F)
    } )
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



################################   Formatting and organizing haplotype frequencies to plot   ###############################
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
hapDTs <- fread("Jan2022.allhaps.restructured.txt.gz", header=TRUE)
lowCovDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Genomewide_avg_coverages/low_cov_samples_DT.txt"), header = T)
lowCovDT[, id := gsub("_", "", id)]
founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)
hapDT <- hapDTs[!lowCovDT, on = .(poolroot = id)]
findNAs <- hapDT[founderfreqs %like% "NA", unique(poolroot)]
hapDT <- hapDT[!poolroot %in% findNAs]
hapsKeepDT <- hapDT[poolroot %in% c("BAS02", "SEE12B02CP010R04", "SEE12B02DI027R01", "SEE12B02YP000R03")]
hapFreqs <- hapsKeepDT[, tstrsplit(founderfreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
names(hapFreqs) <- founderNames
hapsKeepDT$cp <- paste0(hapsKeepDT$chr, "_", hapsKeepDT$pos)
setnames(hapsKeepDT, old = "poolroot", new = "treatment")
metaData <- hapsKeepDT[, c("id", "gp", "chr", "pos", "treatment", "cp")]
allData <- cbind(metaData, hapFreqs)

################################   Formatting and organizing snp frequencies and coverages to plot   ###############################
``
snpDT <- fread("SNPtable.Dec28.sort_ind_restructured.txt", header = T)
snpDF <- as.data.frame(snpDT, stringsAsFactors = FALSE)
snpCovDT <- fread("SNPtable.Dec28.sort_coverage_restructured.txt", header = T)
snpCovDF  <- as.data.frame(snpCovDT, stringsAsFactors = FALSE)

nas <- snpDF[rowSums(is.na(snpDF[,4:21])) > 0,]
nas_idx <- as.numeric(rownames(nas))
if(length(nas_idx) > 0) {
    snpDF2 <- snpDF[-nas_idx,]
    snpCovDF2 <- snpCovDF[-nas_idx,]
    } else {
        snpDF2 <- snpDF
        snpCovDF2 <- snpCovDF
    }
privateSnpDF <- read.delim('Private_snp_lookup_table.txt', header = T)
privateSnpDF <- remove_factors(privateSnpDF)

################################   Making tables of haplotype frequencies vs snp frequencies above a certain coverage cutoff   ###############################

treatments <- c("^BAS02$", "^SEE12B02CP010R04", "^SEE12B02DI027R01", "^SEE12B02YP000R03")
merged_founders <- c('B6')
collapsed_founders <- c('AB3')
bottom_cov_cutoff <- 100
#top_cov_cutoff <- 3000

samples2 <- split(allData, allData$treatment)
snp_dfs <- repeating_list(snpDF2, length(samples2))
snp_cov_dfs <- repeating_list(snpCovDF2 , length(samples2))
private_snp_dfs <- repeating_list(privateSnpDF, length(samples2))
merged_founders <- repeating_list(merged_founders, length(samples2))
collapsed_founders <- repeating_list(collapsed_founders, length(samples2))
bottom_cov_cutoff <- repeating_list(bottom_cov_cutoff, length(samples2))
#top_cov_cutoff <- repeating_list(top_cov_cutoff, length(samples))

hap_calls_vs_snp_calls_bas02_downsampled_df <- Map(private_snp_freqs_vs_hap_freqs, treatments, samples2, snp_dfs, snp_cov_dfs, private_snp_dfs, merged_founders, collapsed_founders, bottom_cov_cutoff)

sample_names <- unlist(lapply(samples2, function(x) {x$treatment[1]}))

saving_lists_of_dfs_to_txt_files(list_of_dfs = hap_calls_vs_snp_calls_bas02_downsampled_df, file_prefix = 'snps_vs_', sample_names = sample_names)

#########################  Analyzing error rate of haplotype caller and down-sampled data  ################
samples <- c("Base population", "Chlorpromazine R04 population wk12", "Diamide R01 population wk12", "YPD R03 population wk12")

files <- dir(pattern = "snps_vs_")
bas02_haps_snps_dfs <- lapply(files, function(x) {
    read.table(x, header = T, sep = "\t")
} )

window <- 10000
error_rate <- 0.05
private_snp_df <- read.delim('Private_snp_lookup_table.txt', header = T)
private_snp_df <- remove_factors(private_snp_df)
windows <- repeating_list(window, length(bas02_haps_snps_dfs))
error_rates <- repeating_list(error_rate, length(bas02_haps_snps_dfs))
private_snp_dfs <- repeating_list(private_snp_df, length(bas02_haps_snps_dfs))

bad_snps <- Map(find_snps_to_filter_out, window = windows, error_rate = error_rates, private_snp_df = private_snp_dfs, bas02_private_snp_df = bas02_haps_snps_dfs)

deleting_snps <- delete_snps(bad_snps)

bas02_filtered_snps <- Map(deleting_snps, bas02_haps_snps_dfs)

counter <- 0
bas02_private_snps <- lapply(bas02_filtered_snps, function(x) {
    counter <<- counter + 1
    x$use <- "no"
    x$shape <- "Unique founders"
    x$use[x$bins ==1] <- "yes"
    x$shape[x$chr == 17] <- "Mitochondrial SNPs"
    x$shape[x$bins == 2] <- "Grouped founders"
    x$sample <- samples[counter]
    x
} )

#snps_vs_haps <- Map(plot_hap_vs_snp_freq, bas02_filtered_snps, bas02_filtered_snps_nogrps)
#########################  Plotting error rate of haplotype caller and down-sampled data  ################

my.formula <- y ~ x

counter <- 0
panelPlots <- lapply(c(bas02_private_snps, bas02_private_snps[length(bas02_private_snps)]), function(snps) {
    counter <<- counter + 1
    snps <- as.data.table(snps)
    snps[, color := shape]
    gplotting <- ggplot(snps[shape == "Unique founders"], aes(unique_snp_freq, hap_calls)) + xlab("SNP frequency") + ylab("Haplotype frequency") + geom_smooth(data = snps[use == "yes"], method = "lm", se=FALSE, colour = "red", size = 1) + stat_poly_eq(data = snps[use == "yes"], formula = my.formula, aes(label = ..rr.label..), parse = TRUE, label.x = 0.5, label.y = 0.75, size = 4, rr.digits = 2) + geom_point(data = snps, aes(unique_snp_freq, hap_calls, colour = as.factor(color)), size = 0.5, alpha = 0.1) + scale_colour_manual(values = c("purple", "black", "gray")) + theme_bw(base_size = 14) + theme(panel.grid = element_blank()) + coord_cartesian(ylim = c(-0.01, 1), xlim = c(-0.01, 1)) + annotate(geom = 'text', label = samples[counter], x = 0, y = Inf, hjust = 0, vjust = 2, size = 4)
    legendPlot <- get_legend(gplotting + guides(colour = guide_legend(ncol = 1, title = "", override.aes = list(size = 3.5, alpha = 1))))
    gplotting2 <- gplotting + theme(legend.position="none")
    if(counter == length(bas02_private_snps)+1) {
        return(legendPlot)
    } else{
        return(gplotting2)
    }
} )

savePlot <- plot_grid(plotlist = panelPlots[-length(panelPlots)], align = 'vh', labels = c("A", "B", "C", "D"), ncol = 2, hjust = -1)
legendAdjust <- plot_grid(panelPlots[[length(panelPlots)]], NULL, ncol = 1)
addLegend <- plot_grid(savePlot, legendAdjust, nrow = 1, rel_widths = c(1,0.2))

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS4_haplotype_caller_validation.pdf"), addLegend, width = 10, height = 10, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS4_haplotype_caller_validation.png"), addLegend, width = 10, height = 10, units = "in", dpi = 350)

########################  Trouble-shooting ######################################