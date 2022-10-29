##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2021
# Title:       SEE01 paper 1 repeatability plots
# Description: Also include scatterplot of LOD scores vs correlation across all positions for all chemicals, fit a lm, then try doing each chemical separately to see if there are some chemicals where there is a correlation b/w LOD score and per-site correlation. May want to try to weight average correlations across replicates so that haplotypes that moved more are given more weight. Maybe can compute mean correlation across replicates for the top two or three moving haplotypes.

##############################################################################

# ============================================================================
# Load packages and sourced files
# Sourced files are kept in the default working directory	

library(tictoc)
library(data.table)
library(tidyverse)
library(scales)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggbeeswarm)
library(DescTools)
library(lemon)
library(abind)
library(job)
library(ggdendro)
library(rtracklayer)
library(GenomicRanges)
library(cowplot)
library(ggpmisc)

source('formatting/Haplotype_file_splitter.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')
source('calculating/Calculating_test_statistics.R')
source('formatting/Haplotype_file_reformatter.R')
source('seq/Position_offsetter.R')



# ============================================================================
# Set global options

defDir <- getwd()
#projectDir <- "/Users/robertlinder/Dropbox/Long_lab/DXQTL03/Primary_experiments/"
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

rep.list <- function(object, repObject) {
    rep(list(object), length(repObject))
}

inflect <- function(avgDiffDT, colNames, threshold = 1){
    avgDiffDF <- as.data.frame(avgDiffDT, stringsAsFactors = FALSE)
    allCols <- lapply(colNames, function(col) {
        up   <- sapply(1:threshold, function(n) c(avgDiffDF[,col][-(seq(n))], rep(NA, n)))
        down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), avgDiffDF[,col][-seq(length(avgDiffDF[,col]), length(avgDiffDF[,col]) - abs(n) + 1)]))
        a <- cbind(avgDiffDF[,col],up,down)
        minMaxList <- list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
        names(minMaxList) <- c(paste0(col, "_minima"), paste0(col, "_maxima"))
        minMaxList
    } ) 
    allCols
}

add.sig.column <- function(avgDiffDT, filteredMinMaxList, LODcutoff) {
    avgDiffDT$significant <- 0
    minOnly <- unlist(filteredMinMaxList[[1]][1])
    maxOnly <- unlist(filteredMinMaxList[[1]][2])
    minRows <- unlist(minOnly)
    maxRows <- unlist(maxOnly)
    avgDiffDT$significant[minRows] <- 2
    avgDiffDT$significant[maxRows] <- 1
    avgDiffDT[LOD < LODcutoff, significant := 0]
    avgDiffDT
}

lod.support.intervals <- function(LODdrop) {
    function(DT) {
        chromSplit <- split(DT, DT$chr)
        LODintervals <- lapply(chromSplit, function(chromDT) {
            #print(chromDT$chr[1])
            #flush.console()
            chromDT[, row := .I]
            chromDT[, LOD2 := LOD - LODdrop]
            rightLod <- chromDT[chromDT, x.row-i.row, on = .(row > row, LOD <= LOD2), mult="first"]
            leftLod <- chromDT[chromDT, x.row-i.row, on = .(row < row, LOD <= LOD2), mult="last"]
            chromDT$rightBound <- rightLod
            chromDT$leftBound <- leftLod
            if(any(is.na(chromDT$rightBound))) {
                chromDT$rightBound[is.na(chromDT$rightBound)] <- nrow(chromDT) - chromDT$row[is.na(chromDT$rightBound)]}
            if(any(is.na(chromDT$leftBound))) {
                chromDT$leftBound[is.na(chromDT$leftBound)] <- 1 - chromDT$row[is.na(chromDT$leftBound)]}
            chromDT
        } )
        chromDTs <- do.call(rbind, LODintervals)
        chromDTs
    }
}

find.lod.dts <- function(DT, peakGP) {
    DT[, row := .I]
    counter <- 0
    findFlanks <- lapply(peakGP, function(peak) {
        counter <<- counter + 1
        peakRow <- DT[gp == peak, row]
        rowsDown <- DT[gp == peak, rightBound]
        rowsUp <- DT[gp == peak, leftBound]
        lodDT <- DT[row <= (peakRow + rowsDown) & row >= (peakRow + rowsUp)][, Idx := counter]
        lodDT
    } )
    flanksDT <- do.call(rbind, findFlanks)
    flanksDT
}

gff_to_granges <- function(gff_file, tags_to_keep) {
    gffRangedData<-import.gff(gff_file)
    myGranges <-as(gffRangedData, "GRanges")
    names(mcols(myGranges))[14] <- "symbol"
    myGranges2 <- myGranges
    mcols(myGranges2) <- NULL
    myGranges2$symbol <- as.character(myGranges$symbol)
    myGranges2$feature <- as.character(myGranges$type)
    myGranges2$name <- as.character(myGranges$Name)
    tags <- as.character(myGranges$type)
    tag_types <- unique(tags)
    tags_keeping <- which(tags %in% tags_to_keep)
    myGranges2 <- myGranges2[tags_keeping]
    tags2 <- myGranges2$feature
    myGranges2$symbol  <- unlist(lapply(1:length(myGranges2), function(x) {
        if(is.na(myGranges2$symbol[x])) {
            return(myGranges2$name[x]) } else {
                return(myGranges2$symbol[x]) }
    } ) )
    myGranges2
}

# peakTable <- Map(peak.table, lodWins, rep.list(myGranges_plot, lodWins), rep.list(3, lodWins))
# lodWin <- lodWins[[1]]
# myGranges2 <- myGranges_plot
# numPeaks <- 3
peak.table <- function(lodWin, myGranges2, numPeaks) {
    tic('Total')
    findTopPeaks <- lodWin[significant == 1][order(-LOD, )][1:numPeaks][, Idx]
    topPeaksDT <- lodWin[Idx %in% findTopPeaks]
    splitIdxs <- split(topPeaksDT, topPeaksDT$Idx)
    genomeCoordinates <- lapply(splitIdxs, function(win) {
        chrName <- as.roman(win$chr[1])
        chr <- paste0("chr", chrName)
        start <- min(win$pos)
        end <- max(win$pos)
        interval <- end-start
        chem <- gsub("18way_", "", win$Chemical[1])
        chemShort <- treatmentAbbrv$abbr[treatmentAbbrv$chemical == chem]
        maxLOD <- win[, round(max(LOD), 0)]
        coordsDF <- data.frame(chrom = chr, start = start, end = end)
        coordsGR <- makeGRangesFromDataFrame(coordsDF)
        overlappingRegion <- subsetByOverlaps(myGranges2, coordsGR)
        genes <- overlappingRegion$symbol
        geneNum <- length(genes)
        geneDT <- data.table(chemical = chemShort, chr = chrName, start_pos = start, end_pos = end, interval_length = interval/1000, genes = I(list(genes)), geneNum = geneNum, LOD = maxLOD)
        setnames(geneDT, "interval_length", "length(kb)")
    } )
    table1 <- do.call(rbind, genomeCoordinates)
    toc()
    table1		
}

extract.peak.hap.changes <- function(peak, DT) {
    peakDT <- peakTables[ID == peak]
    chr <- peakDT$chr
    startPos <- peakDT$start_pos
    endPos <- peakDT$end_pos
    #chrom <- as.numeric(as.roman(strsplit(chr, "chr")[[1]][2]))
    chrom <- as.numeric(as.roman(chr))
    chemical <- peakDT$chemical
    chemDT <- DT[grepl(chemical, id)][chr == chrom & pos >= startPos & pos <= endPos]
    freqDifs <- chemDT[, tstrsplit(avgDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
    avgFreqs <- as.data.frame(chemDT[, tstrsplit(avgFreq, split = ";", type.convert = TRUE, fixed = TRUE)])
    baseFreq <- as.data.frame(chemDT[, tstrsplit(baseFreqs, split = ";", type.convert = TRUE, fixed = TRUE)])
    collHaps <- as.data.frame(chemDT[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
    maxVals <- apply(freqDifs, 1, max, na.rm = T) ## only looking at increasing values
    maxIdx <- apply(freqDifs, 1, which.max)
    maxHap <- collHaps[cbind(seq_along(maxIdx), maxIdx)]
    maxEvFreq <- avgFreqs[cbind(seq_along(maxIdx), maxIdx)]
    maxBaseFreq <- baseFreq[cbind(seq_along(maxIdx), maxIdx)]
    chemDT[, c("maxHap", "maxChange", "maxEvFreq", "maxBaseFreq") := .(maxHap, maxVals, maxEvFreq, maxBaseFreq)]
    nxtMax <- apply(freqDifs, 1, function(x) rev(sort(x))[2])
    nxtMaxIdx <- apply(freqDifs, 1, function(x) which(x == rev(sort(x))[2])[1])
    maxHap2 <- collHaps[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    maxEvFreq2 <- avgFreqs[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    maxBaseFreq2 <- baseFreq[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    chemDT[, c("maxHap2", "maxChange2", "maxEvFreq2", "maxBaseFreq2") := .(maxHap2, nxtMax, maxEvFreq2, maxBaseFreq2)]
    idx <- chemDT[, c("chr", "pos", "gp", "Chemical", "maxHap", "maxChange", "maxEvFreq", "maxBaseFreq", "maxHap2", "maxChange2", "maxEvFreq2", "maxBaseFreq2")]
    idx[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
    idxMelt1 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("maxChange", "maxChange2"))
    idxMelt2 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("maxHap", "maxHap2"))
    idxMelt3 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("maxEvFreq", "maxEvFreq2"), variable.name = "evType", value.name = "evFreq")
    idxMelt4 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical"), measure.vars = c("maxBaseFreq", "maxBaseFreq2"), variable.name = "baseType", value.name = "baseFreq")
    setnames(idxMelt2, old = c("variable", "value"), new = c("hapType", "haplotype"))
    idxMrge <- idxMelt1[, c("hapType", "haplotype") := .(idxMelt2$hapType, idxMelt2$haplotype)]
    setnames(idxMrge, old = c("variable", "value"), new = c("changeType", "frequencyChangeNorm"))
    idxMrge2 <- idxMrge[, c("evFreq", "baseFreq") := .(idxMelt3$evFreq, idxMelt4$baseFreq)]
    idxMrge2[, ID := peakDT$ID]
    idxMrge2
}

extract.peak.all.hap.base.freqs <- function(peak, DT) {
    peakDT <- peakTables[ID == peak]
    chr <- peakDT$chr
    startPos <- peakDT$start_pos
    endPos <- peakDT$end_pos
    chrom <- as.numeric(as.roman(strsplit(chr, "chr")[[1]][2]))
    chemical <- peakDT$chemical
    chemDT <- DT[grepl(chemical, Chemical)][chr == chrom & pos >= startPos & pos <= endPos]
    avgFreqs <- unlist(chemDT[, tstrsplit(avgFreq, split = ";", type.convert = TRUE, fixed = TRUE)])
    baseFreq <- chemDT[, tstrsplit(baseFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
    collHaps <- chemDT[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
    hapDT <- data.table(chr = unique(chemDT$chr), pos = rep(chemDT$pos, length(avgFreqs)/chemDT[,.N]), gp = rep(chemDT$gp, length(avgFreqs)/chemDT[,.N]), Chemical = unique(chemDT$Chemical), avgFreqs = avgFreqs, baseFreq = unlist(baseFreq), haps = unlist(collHaps))
    hapDT <- na.omit(hapDT)
    hapDT
}

# ============================================================================
# Load data needed for downstream analyses

founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)
indLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Hap_adjusted_LOD_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^SEE12B02")
avgHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Averaged_and_sd_tx_tables/", analysisType = "hap_freqs_avg_sd_difs", samplePattern = ".*")


treatKeyDT <- fread("treatment_key.txt", header = T)
plottingFactors <- read.table("all_chems_list.txt", header = T, sep = "\t")
plottingFactors2 <- plottingFactors$Chemical
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
colPal <- c("dodgerblue2", "blue2", "skyblue3", "darkturquoise",
            "brown", "olivedrab", "tomato",
            "chocolate", "burlywood4",
            "saddlebrown",
            "magenta4", "maroon2",
            "mediumpurple",
            "darkcyan",
            "cornsilk4",
            "mediumvioletred")
plottingColors <- colPal
colorDF <- data.frame(factors = plottingFactors2, colors = plottingColors, stringsAsFactors = FALSE)
myColors <- colorDF$colors
names(colPal) <- plottingFactors2

### Creating dictionary of full chemical names to abbreviations

treatmentAbbrv <- data.table(chemical2 = unique(names(colPal)), chemical = c("dmso", "caffeine", "cadmium_chloride", "chlorpromazine", "cisplatin", "tunicamycin", "fluconazole", "ethanol", "urea", "glacial_acetic_acid", "YPD", "nicotinamide", "nicotine", "diamide", "sodium_chloride", "sodium_sulfite"),abbr = c("DM", "CA", "CD", "CP", "CI", "TU", "FL", "ET", "UR", "GA", "YP", "NM", "NT", "DI", "NC", "SO"))
treatmentAbbrv[, chemical3 := paste0("18way_", chemical)]
treatmentAbbrv[, Chemical := gsub("_", " ", chemical)]

fwrite(treatmentAbbrv, "chemAbbrev.txt", sep = "\t", col.names = T)

treatmentAbbrv <- fread("chemAbbrev.txt")
names(colPal) <- treatmentAbbrv$abbr

# ============================================================================
# Find the most significant peaks, the middle peak, and a smaller peak for each drug. Find the 2.5 LOD support interval surrounding the peak LOD scores to subset these positions for plotting.

gff_file <- "saccharomyces_cerevisiae_R64-2-1_20150113.gff"

tags_to_keep <- c("X_element", "telomeric_repeat", "gene", "ARS", "long_terminal_repeat", "ncRNA_gene", "tRNA_gene", "snoRNA_gene", "centromere", "LTR_retrotransposon", "transposable_element_gene", "pseudogene", "Y_prime_element", "telomerase_RNA_gene", "snRNA_gene", "silent_mating_type_cassette_array", "mating_type_region", "intein_encoding_region", "rRNA_gene", "external_transcribed_spacer_region", "internal_transcribed_spacer_region", "non_transcribed_region", "origin_of_replication")

myGranges_plot <- gff_to_granges(gff_file, tags_to_keep )

colNames <- rep(list('LOD'), length(indLODDTs))
threshold <- rep(list(50), length(indLODDTs))
allPeaks <- Map(inflect, indLODDTs, colNames, threshold)
peakCol <- Map(add.sig.column, indLODDTs, allPeaks, rep.list(5, indLODDTs))
top <- lapply(peakCol, function(x) {x[significant == 1][, Idx := 1:.N]})
peakPos <- lapply(top, function(x) {x$gp})
lodInts <- lod.support.intervals(2.5)
addLodInts <- Map(lodInts, peakCol)
lodWins <- Map(find.lod.dts, addLodInts, peakPos)
peakTable <- Map(peak.table, lodWins, rep.list(myGranges_plot, lodWins), rep.list(3, lodWins))
peakTables <- do.call(rbind, peakTable)
peakTables[, ID := c(paste0(rep("CD", 3), 1:3), paste0(rep("CP", 3), 1:3), paste0(rep("DI", 3), 1:3), paste0(rep("GA", 3), 1:3), paste0(rep("NC", 3), 1:3), paste0(rep("UR", 3), 1:3), paste0(rep("YP", 3), 1:3))]
peakTables2 <- peakTables[order(chr, start_pos)]
peakTables3 <- peakTables2[, c("chemical", "ID", "chr", "start_pos", "length(kb)", "geneNum", "LOD", "genes")]
fwrite(peakTables2, paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/peakTable_v5.txt"), sep = "\t", col.names = T)  

allGenes <- unlist(peakTables$genes)
fwrite(data.table(Genes = allGenes), paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/genes_in_peaks.txt"), sep = "\t", col.names = T)  


avgIntLen <- mean(peakTables$`length(kb)`)
geneNums <- dcast(peakTables[, unlist(genes), by = seq_len(nrow(peakTables))][
    , .N, by = .(seq_len, V1)], 
    seq_len ~ V1, value.var = "N", fill = 0L)
geneTotals <- rowSums(geneNums[, -1])
avgGeneNums <- mean(geneTotals)

hapFreqsDT <- do.call(rbind, avgHapDifsDTs)

peakLst <- peakTables[, unique(ID)]
DTs <- rep.list(hapFreqsDT, peakLst)
findHaps <- Map(extract.peak.hap.changes, peakLst, DTs)
allHapPeaks <- do.call(rbind, findHaps)


allFreqs <- Map(extract.peak.all.hap.base.freqs, peakLst, DTs)
allFreqs2 <- do.call(rbind, allFreqs)
#numHaps <- allHapPeaks[, freqBins := cut(avgFreq, seq(0, 1, 0.1))][, freqBins := factor(freqBins, levels = rev(levels(freqBins)))]

# ============================================================================
# Plot the frequency of the average haplotype frequency change of the most changed haplotype vs the next-most changed haplotype 


chemDT <- dcast(allHapPeaks, chr + pos + gp + Chemical + ID ~ changeType, value.var = "frequencyChangeNorm")
dummy <- chemDT[, .(maxChange = max(maxChange)), by = Chemical][, maxChange2 := maxChange]

avgMaxChngeGW <- chemDT[, mean(maxChange)]
avgMaxChnge2GW <- chemDT[, mean(maxChange2)]

maxToNxtMax <- avgMaxChngeGW/avgMaxChnge2GW

chemDTMrge <- merge(chemDT, treatmentAbbrv, by = "Chemical")
chemDT <- chemDTMrge[, Chemical := abbr]
colorPal <- colPal[names(colPal) %in% chemDT$Chemical]
colorPal <- colorPal[order(names(colorPal))]

gplot <- ggplot(data=chemDT, aes(maxChange, maxChange2)) + xlab("Average change of most increased haplotype") + ylab("Average change of 2nd most increased haplotype") + geom_abline(intercept = 0, slope = 1, colour = "black", lty = 3) + geom_smooth(span = 1, colour = "red", size = 0.5, se = FALSE)  + geom_point(aes(colour = Chemical), size = 2, alpha = 0.75) + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), aspect.ratio = 1) + scale_colour_manual(values = colorPal) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure5_1st_2nd_increased_haps_peaks_v2.pdf"), gplot , width = 8, height = 9, units = "in")
system2(command = "pdfcrop", args = c(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure5_1st_2nd_increased_haps_peaks_v2.pdf"), paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure5_1st_2nd_increased_haps_peaks_cropped_v2.pdf")))
# ============================================================================
# Plot the frequency in the base vs the average over replicates of the most changed haplotype and the next-most changed haplotype (different symbols for each) and the frequency in the base vs the average increase in frequency of the most changed and next-most changed haplotypes

allHapPeaks[Chemical == "YPD", Chemical := "ypd"]
#colPalUsing <- colPal[names(colPal) %in% allHapPeaks$Chemical]

colPalUsing  <- colPal[names(colPal) %in% allHapPeaks$Chemical]
colPalUsing  <- colorPal[order(names(colorPal))]
names(colPalUsing) <- unique(allHapPeaks$Chemical)

my.formula <- y ~ x

gplot2 <- ggplot(data=allHapPeaks[changeType == "maxChange"], aes(baseFreq, frequencyChangeNorm)) + xlab("Base frequency") + ylab("Average frequency change") + geom_point(aes(colour = Chemical), size = 2, alpha = 0.75) + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), aspect.ratio = 1) + scale_colour_manual(values = colPalUsing) + coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 1)) + geom_smooth(method = "lm", se=FALSE, colour = "red", size = 0.5) + stat_poly_eq(formula = my.formula, aes(label = ..rr.label..), parse = TRUE, label.x = 1, label.y = 1, size = 3, rr.digits = 2)
legendPlot2 <- get_legend(gplot2 + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))))
gplot2 <- gplot2 + theme(legend.position="none")

gplot3 <- ggplot(data=allHapPeaks[changeType == "maxChange2"], aes(baseFreq, frequencyChangeNorm)) + xlab("Base frequency") + ylab("Average frequency change") + geom_point(aes(colour = Chemical), size = 2, alpha = 0.75) + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position="none") + scale_colour_manual(values = colPalUsing) + coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 1)) + geom_smooth(method = "lm", se=FALSE, colour = "red", size = 0.5) + stat_poly_eq(formula = my.formula, aes(label = ..rr.label..), parse = TRUE, label.x = 1, label.y = 1, size = 3, rr.digits = 2)

savePlot <- plot_grid(gplot2, gplot3, align = 'vh', labels = c("A", "B"), ncol = 1, hjust = -1)
legendAdjust <- plot_grid(legendPlot2, NULL, ncol = 1)
addLegend <- plot_grid(savePlot, legendAdjust, nrow = 1, rel_widths = c(1,0.35))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS15_1st_2nd_increased_haps_peaks_basefreq_vs_change.pdf"), addLegend, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS15_1st_2nd_increased_haps_peaks_basefreq_vs_change.png"), addLegend, width = 8, height = 9, units = "in", dpi = 350) 

# ============================================================================
# Trouble-shooting