##############################################################################

# Project:     DXQTL02
# Author:      Robert Linder
# Date:        2019
# Title:       title
# Description: A brief description

##############################################################################

# ============================================================================
# Load packages and sourced files
# Sourced files are kept in the default working directory	
library(rtracklayer)
library(Gviz)
library(R.utils)
library(GenomicRanges)
library(tictoc)
library(data.table)
library(tidyverse)
library(scales)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggbeeswarm)
library(ggExtra)
library(DEGreport)
library(DescTools)
library(GISTools)
library(pcaPP)
library(grid)
library(gridExtra)
library(cowplot)
library(gtable)
library(NbClust)
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

makePairs <- function(data) 
{
    grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
    #grid <- subset(grid, x != y)
    all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
        xcol <- grid[i, "x"]
        ycol <- grid[i, "y"]
        data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
                   x = data[, xcol], y = data[, ycol], data)
    }))
    all$xvar <- factor(all$xvar, levels = names(data))
    all$yvar <- factor(all$yvar, levels = names(data))
    densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
        data.frame(xvar = names(data)[i], yvar = names(data)[i], x = data[, i])
    }))
    list(all=all, densities=densities)
}

makeUniquePairs <- function(data) 
{
    grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
    grid <- subset(grid, x != y)
    gridDT <- as.data.table(grid)
    dt2 <- unique(gridDT, by=c("x", "y"))[x > y, c("x", "y") := .(y, x)]
    res2 <- unique(dt2, by=c("x", "y"))
    res3 <- as.data.frame(res2)
    all <- do.call("rbind", lapply(1:nrow(res3), function(i) {
        xcol <- res3[i, "x"]
        ycol <- res3[i, "y"]
        data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
                   x = data[, xcol], y = data[, ycol], data)
    }))
    all$xvar <- factor(all$xvar, levels = names(data))
    all$yvar <- factor(all$yvar, levels = names(data))
    densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
        data.frame(xvar = names(data)[i], yvar = names(data)[i], x = data[, i])
    }))
    list(all=all, densities=densities)
}

allUniquePairs <- function(data) 
{
    grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
    grid <- subset(grid, x != y)
    gridDT <- as.data.table(grid)
    dt2 <- unique(gridDT, by=c("x", "y"))[x > y, c("x", "y") := .(y, x)]
    res2 <- unique(dt2, by=c("x", "y"))
    colCombos <- data.frame(apply(res2, 1, function(x) names(data)[x]))
    colCombos
}

swr = function(string, nwrap=16) {
    paste(strwrap(string, width=nwrap), collapse="\n")
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

# ============================================================================
# Load data needed for downstream analyses

privateSNPDT <- fread('Private_snps_AB3B6_lookup_table.txt', header = T)
privateSNPDT$unique_founder[grep("_", privateSNPDT$unique_founder)] <- "AB3"
avgSnpDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Avg_snps_abs_diffs_dfs/", analysisType = "snps_abs_diffs", samplePattern = ".*")
avgSnpDifsDTs2 <- lapply(avgSnpDifsDTs, function(x) {
    x$founder <- NA
    x$founder[x$gp %in% privateSNPDT$gp] <- privateSNPDT$unique_founder
    names(x)[5] <- "control.difs"
    x
})

indLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Hap_adjusted_LOD_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^SEE12B02")
repLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Repeatability_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^18way")
avgHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Averaged_and_sd_tx_tables/", analysisType = "hap_freqs_avg_sd_difs", samplePattern = ".*")

treatKeyDT <- fread("treatment_key.txt", header = T)
topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))
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
gray <- col2rgb('grey50')
chromBarCol <- rgb(gray[1],gray[2],gray[3], maxColorValue = 255, alpha = 100)
colorCodes = c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
hx = sapply(strsplit(colorCodes, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
hx2 = hx[-c(5,21,24)]   
mycols = c(hx2,"#000000")
mypal <- RColorBrewer::brewer.pal(4, "Dark2")
myColPal <- add.alpha(mycols, 0.5)

# ============================================================================
# Plot super-imposed  zscores for chlorp, nacl,and YPD based on transformed LOD scores

chroms <- as.character(as.roman(1:16))
popDT <- do.call(rbind, indLODDTs)
popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
lodDT <- dcast(popDT, gp + chr + pos ~ Chemical, value.var = "logLOD")
names(lodDT) <- gsub("18way_", "", names(lodDT))
names(lodDT) <- gsub("_", " ", names(lodDT))
supDT <- lodDT[, c("gp", "chr", "pos", "chlorpromazine", "sodium chloride", "YPD"), with = FALSE]
lodDF <- as.data.frame(supDT)
lodDF[,4:ncol(lodDF)] <- scale(lodDF[,4:ncol(lodDF)])
newDT <- as.data.table(lodDF)

replodDT2 <- melt(newDT, id.vars = c("gp", "chr", "pos"), measure.vars = c("chlorpromazine", "sodium chloride", "YPD"),variable.name = "group", value.name = "zscore")
replodDT2[, chrBounds := "gray"][(chr %% 2) == 0, chrBounds := "white"]

start <- ch.bounds[seq(1, length(ch.bounds)-2, 2)]
end <- ch.bounds[seq(2, length(ch.bounds)-2, 2)]
chrDT <- data.table(start = start, end = end, chrBounds = rep("gray", 8))

fig7 <- ggplot() + geom_line(data = replodDT2, aes(gp, zscore, colour = group, group = group)) + geom_rect(data = chrDT, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.3) + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), legend.position = "bottom") + scale_colour_manual(values = alpha(c("blue", "red", "black"), 0.4)) + scale_x_continuous(limits = c(0, ch.bounds[17]), breaks = c(g_l[1:16] + offsets[[2]][1:16]/2), labels = as.roman(1:16), expand = c(0, 0)) + coord_cartesian(ylim = c(-5,5)) + xlab("")
fig7updated <- lemon::reposition_legend(fig7, 'top right',  offset=0.002)

# ============================================================================
# Calculating z-scores for log transformed LOD data and extracting rows where multiple chemicals are +/- 1.96 sd from their mean

popDT <- do.call(rbind, indLODDTs)
popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
chroms <- as.character(as.roman(1:16))
filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Pleiotropy_plots/")
lodDT <- dcast(popDT, gp + chr + pos ~ Chemical, value.var = "logLOD")
names(lodDT) <- gsub("18way_", "", names(lodDT))
names(lodDT) <- gsub("_", " ", names(lodDT))
lodDF <- as.data.frame(lodDT)
lodDF[,4:ncol(lodDF)] <- scale(lodDF[,4:ncol(lodDF)])
chemDF <- lodDF[,4:ncol(lodDF)]
chemDF[chemDF < 1.96] <- 0
#chemDF[chemDF < 1.65] <- 0
chemDF[chemDF != 0] <- 1
chemDF$sigRow <- rowSums(chemDF)
lodDF$sigRow <- chemDF$sigRow
#lodDF[,4:ncol(lodDF)] <- abs(lodDF[,4:ncol(lodDF)])
supDT <- as.data.table(lodDF)
chems <- names(supDT)[4:(ncol(supDT)-1)]
supDTcpy <- copy(supDT)
supDT <- supDT[sigRow > 1] ### AT least two drugs share a significantly increased locus
supDF <- as.data.frame(supDT)
frame <- as.matrix(supDF[, 4:(ncol(supDF)-2)])
sigPos <- as.matrix(which(frame>=1.96, arr.ind = TRUE))
sigPos[,2] <- sigPos[,2] + 3
sigFreqs <- supDF[sigPos]
sigDT <- data.table(chr = supDT$chr[sigPos[,1]], pos = supDT$pos[sigPos[,1]], gp = supDT$gp[sigPos[,1]], Chemical = names(supDT)[sigPos[,2]], zscore = sigFreqs)[order(gp, Chemical)]
hcl=hclust(dist(sigDT$gp))
#plot(hcl,labels = FALSE, hang= -1)
#rect.hclust(hcl, h = 200000, border = "red")
clu.h100=cutree(hcl,h=200000)
sigDT[, Idx := clu.h100][Idx %in% 2:6, Idx := Idx+1][Idx == 1 & pos > 595000, Idx := 2]
sigDTLst <- split(sigDT, sigDT$Idx)
maxInt <- lapply(sigDTLst, function(x) {
    gpRange <- unique(supDTcpy[chr == x$chr[1] & gp >= (x$gp[1] - 5000) & gp <= (x$gp[length(x$gp)] + 5000), gp])
    chemNum <- sigDT[Idx == x$Idx[1], uniqueN(Chemical)]
    chems <- sigDT[Idx == x$Idx[1], unique(Chemical)]
    rangeDT <- data.table(gp = rep(gpRange, chemNum), chr = x$chr[1], Chemical = sort(rep(chems, length(gpRange))), Idx = x$Idx[1])
    rangeDT
})
allInts <- do.call(rbind, maxInt)
zscoreDT <- melt(supDTcpy, id.vars = c("gp", "chr", "pos"), measure.vars = c(4:(ncol(supDTcpy)-1)), variable.name = "Chemical", value.name = "zscore")
mrgeInts <- merge(allInts, zscoreDT, by = c("gp", "chr", "Chemical"))[, "pos(kb)" := pos/1000][, "pos" := NULL]

# ============================================================================
# Data wrangle the average frequency change of all haplotypes (for each chemical), with each chemical a different shape and the haplotypes colored consistently with previous plots.

chemShapes <- data.table(chems = unique(popDT$Chemical), shapes = c(0:6))
chemShapes2 <- chemShapes[, chems := gsub("18way_", "", chems)][, chems := gsub("_", " ", chems)]

avgDifsDT <- do.call(rbind, avgHapDifsDTs)
avgDifsDT[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
mergeAvgDifs <- merge(avgDifsDT, allInts, by = c("gp", "chr", "Chemical"))
idxAvgSplit <- split(mergeAvgDifs, mergeAvgDifs$Idx)

idxLooperAll <- lapply(idxAvgSplit, function(idx)  {
    chemSplit <- split(idx, idx$Chemical)
    chemLooper <- lapply(chemSplit, function(chem) {
        freqDifs <- chem[, tstrsplit(avgDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
        collHaps <- chem[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
        allHaps <- cbind(chem$chr, chem$pos, chem$gp, chem$Chemical, chem$Idx, collHaps, freqDifs)
        freqsDF <- as.data.frame(freqDifs)
        newFreqs <- lapply(freqsDF, function(y) gsub("^-\\d.*|^\\d.*", 1, y))
        newFreqsDF <- as.data.frame(do.call(cbind, newFreqs))
        newFreqsDF2 <- sapply(newFreqsDF, as.numeric)
        allHaps$reps <- rowSums(newFreqsDF2, na.rm = T)
        positions <- rep(unlist(allHaps[,2]), allHaps$reps) 
        gpositions <- rep(unlist(allHaps[,3]), allHaps$reps)
        allfreqs <- as.numeric(unlist(t(freqDifs)))
        allfreqs2 <- na.omit(allfreqs)
        hapVec <- as.character(unlist(t(collHaps)))
        hapVec2 <- na.omit(hapVec)
        idxDT <- data.table(chr = allHaps$V1[1], pos = positions, gp = gpositions, Chemical = allHaps$V4[1], Idx = allHaps$V5[1], haps = hapVec2, freqDifs = allfreqs2)
        hapColors <- match(idxDT$haps, topHaps)
        hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
        Colors <- mycols[hapColors]
        chemicals <- match(idxDT$Chemical, chemShapes2$chems)
        Shapes <- chemShapes2$shapes[chemicals]
        idxDT[, c("color.codes", "Shape") := .(Colors, Shapes)][, "pos(kb)" := pos/1000][, "pos" := NULL]
        chem <- idxDT
        chem <- chem[, hapGrps := haps][!hapGrps %in% topHaps, hapGrps := "other"][order(haps, gp)][, c("difCol") := .(c(1000, diff(gp))), by = "haps"][, row := .I]
        chem[, grpCol := paste0(haps, "_", difCol)]
        chem2 <- chem[, grp := cumsum(c(TRUE, diff(gp)!=1000)), by = haps][order(row)] ### troubleshoot have two different haps in same group (gp 4207890)
        chem2[, grp := cumsum(c(TRUE, diff(gp)!=1000)), by = haps]
        chem2[, grp2 := rleid(haps)]
        chem2[, grpCor := c(0, as.numeric((diff(grp) != 0)))]
        chem2[, grpCor2 := c(0, as.numeric((diff(grp2) != 0)))]
        chem2[, grpSum := rowSums(chem2[, c("grpCor", "grpCor2")])]
        #rleGrpSum <- rle(chem2$grpSum)
        
        rowCounter <- 0
        grpCounter <- 1
        makeGrps <- lapply(chem2$grpSum, function(x) {
            rowCounter <<- rowCounter + 1
            if(x != 0) {
                grpCounter <<- grpCounter + 1
                return(data.table(row = chem2$row[rowCounter], grp = grpCounter))}
        } )
        grpDT <- do.call(rbind, makeGrps)
        firstRow <- data.table(row = 1, grp = 1)
        lastRows <- data.table(rows = grpDT$row[length(grpDT$row)]:chem2$row[length(chem2$row)], grps = grpDT$grp[length(grpDT$grp)])
        grpMrgeDT <- rbind(firstRow, grpDT)
        expandedDT <- lapply(1:(nrow(grpMrgeDT)-1), function(x) {
            rowIdces <- grpMrgeDT$row[x]:(grpMrgeDT$row[x+1]-1)
            grpFill <- rep(grpMrgeDT$grp[x], length(rowIdces))
            newDT <- data.table(rows = rowIdces, grps = grpFill)
            newDT
        } )
        allGrpDT <- do.call(rbind, expandedDT)
        allGrpMrgeDT <- rbind(allGrpDT, lastRows)
        chem2[, finalGrp := allGrpMrgeDT$grps]
        chem2
    } )
    idxDTs <- do.call(rbind, chemLooper)
    idxGrps <- rle(idxDTs$grp)
    idxGrps$values <- 1:length(idxGrps$values)
    idxDTs$grp <- rep(idxGrps$values, idxGrps$lengths)
    idxDTs
} )    

# ============================================================================
# plot LOD over avg haplotype change for manually specified windows above

idx <- idxLooperAll[[6]]
zDT <- mrgeInts[Idx == idx$Idx[1]]
chemicals <- match(zDT$Chemical, chemShapes2$chems)
shapes <- chemShapes2$shapes[chemicals]
zDT[, Shape := shapes]
posMin <- 165
posMax <- 180
zDT2 <- zDT[`pos(kb)` >= posMin & `pos(kb)` <= posMax]
idx2 <- idx[`pos(kb)` >= posMin & `pos(kb)` <= posMax]
#ylims <- c(0, max(lodDT$LOD))
zPlot <-  ggplot(data=zDT2, aes(`pos(kb)`, zscore, shape = Chemical, group = Chemical)) + geom_point(size = 2) + geom_line() + scale_shape_manual(values = setNames(zDT2$Shape, zDT2$Chemical)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", plot.margin = margin(b = 0, r = 1)) + ylab("zscore") + xlab("") 
avgHapPlot <- ggplot(data=idx2, aes(`pos(kb)`, freqDifs, colour = hapGrps, shape = Chemical, group = interaction(finalGrp, Chemical))) + coord_cartesian(ylim = c(-1, 1)) + geom_point(size = 2) + geom_line() + scale_colour_manual(values=setNames(idx2$color.codes, idx2$hapGrps)) + scale_shape_manual(values = setNames(idx2$Shape, idx2$Chemical)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), plot.margin = margin(t = 0, r = 1), legend.position = "bottom") + ylab("mean haplotype change") + xlab(paste0("chr", as.roman(idx2$chr[1]), " position (kb)"))
legendPlot <- get_legend(avgHapPlot + guides(color = guide_legend(nrow = 2)) + theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin()))
avgHapPlot <- avgHapPlot + theme(legend.position="none")

# ============================================================================
# Plot shared significant peaks gviz; first edit the mrgeInts variable to keep only indices/positions determined by manually exploring the data above

idx7 <- mrgeInts[Idx == 6 & `pos(kb)` >= 165 & `pos(kb)` <= 180][, Idx := 7]
gff_file <- "saccharomyces_cerevisiae_R64-2-1_20150113.gff"

tags_to_keep <- c("X_element", "telomeric_repeat", "gene", "ARS", "long_terminal_repeat", "ncRNA_gene", "tRNA_gene", "snoRNA_gene", "centromere", "LTR_retrotransposon", "transposable_element_gene", "pseudogene", "Y_prime_element", "telomerase_RNA_gene", "snRNA_gene", "silent_mating_type_cassette_array", "mating_type_region", "intein_encoding_region", "rRNA_gene", "external_transcribed_spacer_region", "internal_transcribed_spacer_region", "non_transcribed_region", "origin_of_replication")

myGranges_plot <- gff_to_granges(gff_file, tags_to_keep )
chemShapes <- data.frame(chems = unique(popDT$Chemical), shapes = c(0:6))

avgDifs <- do.call(rbind, avgHapDifsDTs)
plotWinsList <- idx7
fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/")
statistic <- "zscore"

x <- plotWinsList
xCst <- dcast(x, chr + `pos(kb)` + gp ~ Chemical, value.var = "zscore")
xCst[, pos := 1000*`pos(kb)`]
if(xCst$chr[1] == 'chrmt') {chr_name <- 'M'} else {
    chr_name <- as.roman(xCst$chr[1])}
chr <- chr_name
start <- min(xCst$pos)
end <- max(xCst$pos)
grtrack <- GeneRegionTrack(myGranges_plot, chromosome=paste0('chr',chr), start=start, end=end, showId=TRUE, shape='fixedArrow', thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5", "miRNA", "lincRNA"), showFeatureId=TRUE, showTitle = TRUE, name = 'Genes', just.group = "below", geneSymbol = TRUE, feature = myGranges_plot$feature, featureAnnotation = unique(myGranges_plot$feature), arrowHeadWidth=10, cex = 0.5)
pdf(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure8C.pdf"), height = 1, width = 8)
genePlot <- plotTracks(list(grtrack), from = start, to = end, cex.axis = 1, col.axis = "black", chromosome = paste0('chr',chr), showTitle = FALSE )
dev.off()

savePlot3 <- plot_grid(fig7updated, zPlot, avgHapPlot, align = 'vh',labels = c("A", "B", "D"), ncol = 1, hjust = -1)
addLegend <- plot_grid(savePlot3, legendPlot, ncol = 1, rel_heights = c(1, .1))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/Figure9ABD.pdf"), addLegend, width = 8, height = 9, units = "in")

# ============================================================================
# Trouble-shooting