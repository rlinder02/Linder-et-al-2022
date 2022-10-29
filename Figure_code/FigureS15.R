##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2019
# Title:       Plotting raw haplotype frequencies 
# Description: Plotting haplotype frequencies of the base population, along with the treatment replicates and ypd-only replicates; also making a table of raw haplotype frequencies and heterozygosity 

##############################################################################

# ============================================================================
# Load packages and sourced files
# Sourced files are kept in the default working directory 
library(gprofiler2)
library(gridExtra)
library(egg)
library(readxl)
library(scales)
library(seqinr)
library(tictoc)
library(data.table)
library(R.utils)
library(dplyr)
library(ggplot2)
library(GGally)
library(gt)
library(ggpubr)
library(ggcorrplot)
library(ggbeeswarm)
library(RColorBrewer)
library(gridGraphics)
library(lemon)
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')
source('seq/Position_offsetter.R')

# ============================================================================
# Set global options

defDir <- getwd()
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

# ============================================================================
# Custom functions

add.alpha <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
              rgb(x[1], x[2], x[3], alpha=alpha))  
}

# ============================================================================
# Load data
indLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Hap_adjusted_LOD_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^SEE12B02")
indHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_het_diffs_tables/", analysisType = "haps_het_diffs", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]_")
treatmentDT <- fread("SEE01_unique_samples.txt", header = F)
founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
treatKeyDT <- fread("treatment_key.txt", header = T)
offsets <- read.table("newoffsets.txt", header = T)
offsets$lines <- ceiling(offsets[,2]/50)
offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
offsets$chr <- 1:17
g_l <- c(0, cumsum(offsets$len))
ch.bounds <- c(0, g_l[1:17] + offsets[,2])
#correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
max.pos <- g_l[length(g_l)]
mid.ch <- diff(ch.bounds)/2
midpt.ch <- ch.bounds[2:18] - mid.ch
plottingFactors <- founderNames
colPal <- c("dodgerblue2", "blue2", "skyblue3", "darkturquoise",
            "brown", "olivedrab", "tomato",
            "chocolate", "burlywood4",
            "saddlebrown",
            "magenta4", "maroon2",
            "mediumpurple",
            "darkcyan",
            "cornsilk4",
            "mediumvioletred",
            "slategray3")
plottingColors <- colPal
colorDF <- data.frame(factors = plottingFactors, colors = plottingColors, stringsAsFactors = FALSE)
myColors <- colorDF$colors
myColorsAlpha <- add.alpha(myColors, alpha = 1)

# ============================================================================
# Plotting mean per-site LOH vs per-site LOD for all chemicals (separate panels for each chemical)

repsUsing <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"))
repsUsing[, chemical := gsub("_12", "", chemical)]
popDT <- do.call(rbind, indHapDifsDTs)
popsUsing <- popDT[Chemical %in% unique(repsUsing$chemical)]
popsUsing[, LOH := baseHeterozygosityCol - heterozygosity]
avgLOH <- popsUsing[, .(avgLOH = mean(LOH)), by = .(Chemical, gp)]

lodDT <- do.call(rbind, indLODDTs)
chem2 <- merge(avgLOH, lodDT, by = c("gp", "Chemical"))
chem2[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
chemDF <- as.data.frame(chem2)
panelPlot <- ggplot(chemDF, aes(LOD, avgLOH)) + facet_wrap(~as.factor(Chemical), scales = "free_x", labeller = labeller(xvar = label_wrap_gen(16))) + geom_hex(bins = 50) + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "right", label.y = 0.75, size = 3, hjust = 1) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3, se=FALSE, colour = "red", size = 0.5) + theme_bw(base_size = 10) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) 
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/FigureS14_LODvsLOH_genomewide.pdf"), panelPlot, width = 8, height = 10, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/FigureS14_LODvsLOH_genomewide.png"), panelPlot, width = 8, height = 10, units = "in", dpi = 350)



# ============================================================================
# Trouble-shooting