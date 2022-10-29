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
treatmentAbbrv <- fread("chemAbbrev.txt")

# ============================================================================
# Plotting summary of mean per-site heterozygosity change for each replicate

repsUsing <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"))
allBind <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_All_Reps_Hets.txt"))
allBind <- allBind[chemical != "18way_base_2"][order(chemical)]
allBind <- allBind[chemical %in% unique(repsUsing$chemical)]
allBind <- allBind[, chemical := gsub("_12|18way_", "", chemical)]
allBind <- allBind[, chemical := treatmentAbbrv$abbr[treatmentAbbrv$chemical == chemical], by = chemical][order(chemical, Replicate)]
repsUsing <- repsUsing[, chemical := gsub("_12|18way_", "", chemical)]
repsUsing <- repsUsing[, chemical := treatmentAbbrv$abbr[treatmentAbbrv$chemical == chemical], by = chemical][order(chemical, Replicate)][, status := "Included"]
allBind[reason == "correlation", reason := "Correlation excluded"]
allBind[reason == "heterozygosity", reason := "Heterozygosity excluded"]
allBind2 <- allBind[!(id %in% repsUsing$id)]
setnames(allBind2, old = "reason", new = "status")
allReps <- rbind(repsUsing, allBind2)

hapDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"), header = T)
hapDT$treatment <- gsub("^18way_|_([0-1]?[0-9])$", "", hapDT$chemical)
hapDT <- hapDT[, treatment := treatmentAbbrv$abbr[treatmentAbbrv$chemical == treatment], by = treatment][order(treatment)]
nReps <- hapDT[, uniqueN(Replicate), by = treatment]
hapDT2 <- merge(hapDT, nReps, by = "treatment")
hapDT2[, treatment := paste0(treatment, "\n", "(", V1, " reps)")]
hetDT <- hapDT[, .(avgHetChange = mean(hetChange)), by = treatment]
baseDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_All_Reps_Hets.txt"), header = T)
baseDT <- baseDT[id == "BAS02"]
avgBaseDT <- baseDT[, .(het = mean(het)), by = id]
avgBaseDT <- cbind(avgBaseDT, baseDT[1, c("chemical", "Replicate", "hetChange")])
#hapDT$treatment <- as.factor(hapDT$treatment)
hetsExc <- allReps[status == "Heterozygosity excluded", .N, by = chemical]
countReps <- allReps[, .N, by = chemical]
fracExc <- hetsExc$N/countReps$N
## Creating a dictionary of full chemical names to short names


### Panel A

repsUsingPlot <- ggplot(allReps, aes(chemical, hetChange, colour = status)) + geom_point(size = 2) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + scale_colour_manual(values = c("orange", "red", "black")) + coord_cartesian(ylim = c(0, 0.25)) + scale_y_continuous(name = "mean per-site heterozygosity deviation", breaks = seq(0, 0.25, 0.05), labels = seq(0, 0.25, 0.05)) + xlab("")
revisedPlot <- reposition_legend(repsUsingPlot, 'top left',  offset=0.002)

### Panel B

hetPlot <- ggplot(hapDT2, aes(x=treatment, y=het)) + geom_boxplot(notch=FALSE) + geom_beeswarm(size = 1, color = "blue") + theme_bw(base_size = 10) + geom_hline(yintercept = avgBaseDT$het, linetype = "dashed", color = "red") + coord_cartesian(ylim = c(0, 0.8)) + xlab("") + ylab("heterozygosity")

fig2 <- ggarrange(revisedPlot, hetPlot, ncol = 1, labels="AUTO", align = "h")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS7_rep_hets.png"), fig2, width = 8.5, height = 9, units = "in", dpi = 350)
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS7_rep_hets.pdf"), fig2, width = 8.5, height = 9, units = "in")

# ============================================================================
# Trouble-shooting