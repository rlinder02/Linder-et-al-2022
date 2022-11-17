##############################################################################

# Project:     DXQTL02
# Author:      Robert Linder
# Date:        2019
# Title:       title
# Description: Calculating the sum of squared differences between the haplotype frequencies of all samples and bas02 as well as average per-site heterozygosity differences; storing the replicates not using in the Ind_tx_ind_haps_sum_of_squared_diffs_tables_deprecated folder.

##############################################################################

# ============================================================================
# Load packages and sourced files
# Sourced files are kept in the default working directory	
#BiocManager::install("preprocessCore")
library(tictoc)
library(data.table)
source('formatting/Haplotype_file_splitter.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')
source('calculating/Calculating_test_statistics.R')
source('formatting/Haplotype_file_reformatter.R')
source('seq/Position_offsetter.R')



# ============================================================================
# Set global options

defDir <- getwd()
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

# ============================================================================
# Custom functions

rep.list <- function(object, repObject) {
    rep(list(object), length(repObject))
}
# ============================================================================
# Load data needed for downstream analyses

founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)
lowCovDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Genomewide_avg_coverages/low_cov_samples_DT.txt"), header = T)
lowCovDT[, id := gsub("_", "", id)]

indHapDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Individual_tx_tables/", analysisType = "hap_freqs", samplePattern = "^SEE12B02|^BAS02")
indSnpData <- fread("SNPtable.April23.sort_ind_restructured.txt", header=TRUE)
names(indSnpData)[names(indSnpData) == 'A6'] <- 'FA06'
metaData <- fread("sexual.rawnames.txt", sep = "\t", header = F)
repsClass <- fread("all_samples_classified.txt", sep = "\t", header = T)

# ============================================================================
# Reformat the snp frequencies file to host (for all samples, excluding low coverage ones)

snpDataInfo <- indSnpData[, c("chr", "POS", "gp", "REF", "ALT")]
snpDataFoundersBase <- indSnpData[, c("AB1b", "AB2b", "AB3b", "AB4b", "A5", "FA06", "A7", "A8", "A9", "A11", "A12", "B5", "B6", "B7", "B8", "B9", "B11", "B12", "BAS02")]
setnames(snpDataFoundersBase, old = c("AB1b", "AB2b", "AB3b", "AB4b", "FA06", "BAS02"), new = c("AB1", "AB2", "AB3", "AB4", "A6", "18F12v2"))
indSnpDataDF <- as.data.frame(indSnpData)
snpDataEvolved <- indSnpDataDF[, grep("SEE12B02", names(indSnpDataDF))]
snpDataEvolvedFilt <- snpDataEvolved[, !names(snpDataEvolved) %in% lowCovDT$id]

snpData <- cbind(snpDataInfo, snpDataFoundersBase, snpDataEvolvedFilt)

fwrite(snpData, "sexual_paper_snp_frequencies.txt", sep = "\t", col.names = T)
# ============================================================================
# Reformat the collapsed haplotype frequency file to host (for all samples that passed low coverage filter)

popDT <- do.call(rbind, indHapDTs)
popDT <- popDT[!lowCovDT, on = "id"]
popDT <- popDT[, c("chr", "pos", "gp", "id", "collapsedFreqs", "collapsedFounders", "heterozygosity")]
popDT[id == "BAS02", id := "18F12v2"]

fwrite(popDT, "sexual_paper_haplotype_frequencies.txt", sep = "\t", col.names = T)
# ============================================================================
# Make a metadata table for all samples, including names of fastq files and status (included, excluded, why excluded); DMSO R10 is right on border of AD vs S (has 5 including this one, plus one that doesn't cluster); may need to add DMSO to figures 

metaData <- metaData[, c("V1", "V3", "V4")]
setnames(metaData, old = c("V1", "V3", "V4"), new = c("id", "file1", "file2"))
metaData2 <- metaData[grepl("SEE12B02", id)][order(id)]
sexualBadClust <- paste0("SEE12B02", c("FLR13", "FLR03", "NMR01", "YPR10", "NTR02", "GAR08", "FLR08", "NTR01", "FLR07", "NTR07", "FLR16", "NMR10", "NMR04", "NMR11", "NMR08", "CDR01", "CDR02", "CDR03", "CPR12", "CPR13", "FLR09", "NMR09", "DMR08"))
chemsExluded <- c("CA", "CI", "DM", "ET", "FL", "NM", "NT", "SO", "TU")
## 3-letter status codes includes xLC (excluded low coverage), xAH (excluded, aneuploid haploid), xCD (excluded, clonal diploid), xPC (sexual, excluded for poor clustering), xFR (sexual, passed all filters, but chemical has too few replicates), xLH (excluded, lower heterozygosity than most other reps), iSP (included, sexual pass; passes all filters)

metaData2[!id %in% popDT$id, Status := "xLC"]
metaData2[, idShort := gsub("[0-9][0-9][0-9]", "", id)]
mrgeData <- merge(metaData2, repsClass[, c("id", "Type")], by.x = "idShort", by.y = "id", all.x = T)
#mrgeData[Status == "xLC", Type := "xLC"]
mrgeData[idShort %in% sexualBadClust & is.na(Status), Status := "xPC"][Status == "xLC", Type := NA]
mrgeData[, chem := substr(id, 9,10)]
mrgeData[chem %in% chemsExluded & Type == "Outbred_sexual" & is.na(Status), Status := "xFR"] ## fix this
mrgeData[is.na(Status) & Type == "Outbred_sexual", Status := "iSP"][Type == "Aneuploid_haploid" & is.na(Status), Status := "xAH"][Type == "Clonal_diploid" & is.na(Status), Status := "xCD"]
mrgeData[id %in% c("SEE12B02NC950R13", "SEE12B02NC950R16", "SEE12B02CP010R05"), Status := "xLH"]
mrgeData[, idShort := NULL]

fwrite(mrgeData, "sexual_paper_sample_metadata.txt", sep = "\t", col.names = T)

## NOTE- changed # outbred sexuals for some (now 105 instead of 110 classified as outbred sexual in general; DMSO still has 6), plus DMSO clustering eliminates one, so only have 5 DMSO reps- so maybe remake dendrogram- but one DMSO rep is right on the edge of heterozygosity- R10); change figures   

# turns out NCR16 and NCR11 and NC13 don't cluster well, but NC16 does;  CPR05 clusters ok with other chlorpromazines 

# ============================================================================
# Trouble-shooting