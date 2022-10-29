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

zscore.support.intervals <- function(zdrop) {
    function(DT) {
        chemSplit <- split(DT, DT$Chemical)
        chemInts <- lapply(chemSplit, function(chem) {
            #print(chromDT$chr[1])
            #flush.console()
            chem[, row := .I]
            chem[, zscore2 := zscore - zdrop]
            rightLod <- chem[chem, x.row-i.row, on = .(row > row, zscore <= zscore2), mult="first"]
            leftLod <- chem[chem, x.row-i.row, on = .(row < row, zscore <= zscore2), mult="last"]
            chem$rightBound <- rightLod
            chem$leftBound <- leftLod
            if(any(is.na(chem$rightBound))) {
                chem$rightBound[is.na(chem$rightBound)] <- nrow(chem) - chem$row[is.na(chem$rightBound)]}
            if(any(is.na(chem$leftBound))) {
                chem$leftBound[is.na(chem$leftBound)] <- 1 - chem$row[is.na(chem$leftBound)]}
            chem
        } )
        chemDTs <- do.call(rbind,  chemInts)
        chemDTs
    }
}

find.lod.dts <- function(DT, peakDT) {
    DT[, row := .I]
    findFlanks <- lapply(1:nrow(peakDT), function(x) {
        counter <- peakDT[x, Idx]
        peakRow <- DT[gp == peakDT[x,gp] & Chemical == peakDT[x,Chemical], row]
        rowsDown <- DT[gp == peakDT[x,gp] & Chemical == peakDT[x,Chemical], rightBound]
        rowsUp <- DT[gp == peakDT[x,gp] & Chemical == peakDT[x,Chemical], leftBound]
        zDT <- DT[row <= (peakRow + rowsDown) & row >= (peakRow + rowsUp)][, Idx := peakDT[x, Idx]]
        zDT[, ID := paste0("PL", counter)]
        zDT
    } )
    flanksDT <- do.call(rbind, findFlanks)
    flanksDT
}

find.idx.range <- function(DT, boundsDT) {
    idxSplit <- split(boundsDT, boundsDT$Idx)
    counter <- 0
    findAll <- lapply(idxSplit, function(idx) {
        counter <<- counter + 1
        chems <- idx[, unique(Chemical)]
        DTint <- DT[gp >= idx$start[1] & gp <= idx$end[1] & Chemical %in% chems][, Idx := idx$Idx[1]][, ID := paste0("PL", counter)]
        DTint
    } )
    flanksDT <- do.call(rbind, findAll)
    flanksDT
}

zscore.sums.support.intervals <- function(zdrop) {
    function(DT) {
        idxSplit <- split(DT, DT$Idx)
        chemInts <- lapply(idxSplit, function(chem) {
            #print(chromDT$chr[1])
            #flush.console()
            setnames(chem, old = "sumZs", new = "zscore")
            chem[, row := .I]
            chem[, zscore2 := zscore - zdrop]
            rightLod <- chem[chem, x.row-i.row, on = .(row > row, zscore <= zscore2), mult="first"]
            leftLod <- chem[chem, x.row-i.row, on = .(row < row, zscore <= zscore2), mult="last"]
            chem$rightBound <- rightLod
            chem$leftBound <- leftLod
            if(any(is.na(chem$rightBound))) {
                chem$rightBound[is.na(chem$rightBound)] <- nrow(chem) - chem$row[is.na(chem$rightBound)]}
            if(any(is.na(chem$leftBound))) {
                chem$leftBound[is.na(chem$leftBound)] <- 1 - chem$row[is.na(chem$leftBound)]}
            chem
        } )
        chemDTs <- do.call(rbind,  chemInts)
        chemDTs
    }
}

find.max.sum.zscores <- function(DT, peaks) {
    DT[, row := .I]
    counter <- 0
    findFlanks <- lapply(peaks, function(peak) {
        peakRow <- DT[gp == peak, row]
        rowsDown <- DT[gp == peak, rightBound]
        rowsUp <- DT[gp == peak, leftBound]
        zDT <- DT[row <= (peakRow + rowsDown) & row >= (peakRow + rowsUp)]
        zDT
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

pleio.peak.table <- function(pleioWins, myGranges2) {
    tic('Total')
    splitIdxs <- split(pleioWins, pleioWins$Idx)
    counter <- 0
    genomeCoordinates <- lapply(splitIdxs, function(win) {
        counter <<- counter + 1
        chrName <- as.roman(win$chr[1])
        chr <- paste0('chr', chrName)
        win[, pos := `pos(kb)`*1000]
        start <- min(win$pos)
        end <- max(win$pos)
        interval <- end-start
        chems <- paste(unique(win$Chemical), collapse = ", ")
        coordsDF <- data.frame(chrom = chr, start = start, end = end)
        coordsGR <- makeGRangesFromDataFrame(coordsDF)
        overlappingRegion <- subsetByOverlaps(myGranges2, coordsGR)
        genes <- overlappingRegion$symbol
        geneDT <- data.table(ID = paste0("PL", counter), chemicals = chems, chr = chr, start_pos = start, end_pos = end, interval_length = interval/1000, genes = I(list(genes)))
        setnames(geneDT, "interval_length", "interval(kb)")
    } )
    table1 <- do.call(rbind, genomeCoordinates)
    toc()
    table1		
}

# ============================================================================
# Load data needed for downstream analyses

indLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Hap_adjusted_LOD_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^SEE12B02")
avgHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Averaged_and_sd_tx_tables/", analysisType = "hap_freqs_avg_sd_difs", samplePattern = ".*")

treatKeyDT <- fread("treatment_key.txt", header = T)
plottingFactors <- read.table("all_chems_list.txt", header = T, sep = "\t")
plottingFactors2 <- plottingFactors$Chemical
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
colorCodes = c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
hx = sapply(strsplit(colorCodes, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
hx2 = hx[-c(5,21,24)]   
mycols = c(hx2,"#000000")
mypal <- RColorBrewer::brewer.pal(4, "Dark2")
myColPal <- add.alpha(mycols, 0.5)
chemShrt <- fread("chemAbbrev.txt")
chemShrt[, chemical2 := gsub("_", " ", chemical)]
fwrite(chemShrt, "chemAbbrev.txt")
# ============================================================================
# Calculating z-scores for log transformed LOD data and extracting rows where multiple chemicals are +/- 1.96 sd from their mean

popDT <- do.call(rbind, indLODDTs)
popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
chroms <- as.character(as.roman(1:16))
filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Pleiotropy_plots/")
lodDT <- dcast(popDT, gp + chr + pos ~ Chemical, value.var = "logLOD")
names(lodDT) <- gsub("18way_", "", names(lodDT))
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
hcl <- hclust(dist(supDT$gp))
plot(hcl,labels = FALSE, hang= -1)
rect.hclust(hcl, h = 150000, border = "red")
clu.h100 <- cutree(hcl,h=150000)
supDT[, Idx := clu.h100]
supDTmlt <- melt(supDT, id.vars = c("gp", "chr", "pos", "Idx"), measure.vars = c("cadmium_chloride", "chlorpromazine", "diamide", "glacial_acetic_acid", "sodium_chloride", "urea", "YPD"), variable.name = "Chemical", value.name = "zscore")
supDTCpymlt <- melt(supDTcpy, id.vars = c("gp", "chr", "pos"), measure.vars = c("cadmium_chloride", "chlorpromazine", "diamide", "glacial_acetic_acid", "sodium_chloride", "urea", "YPD"), variable.name = "Chemical", value.name = "zscore")[, "pos(kb)" := pos/1000][, "pos" := NULL]
sigDT <- supDTmlt[zscore >= 1.96]
chemNums <- sigDT[, .(numChems = uniqueN(Chemical)), by = Idx]
numMrge <- merge(sigDT, chemNums, by = "Idx")
intervals <- numMrge[, .(startPos = rep(min(pos), numChems[1]), stopPos = rep(max(pos), numChems[1]), chr = rep(min(chr), numChems[1]), Chemical = unique(Chemical)), by = Idx]
gpRange <- intervals[, .(pos = startPos:stopPos), by = .(Idx, Chemical, chr)]
rangeZs <- merge(gpRange, supDTmlt, by = c("Idx", "Chemical", "chr", "pos"))[, "pos(kb)" := pos/1000][, "pos" := NULL]
rangeZs$Chemical <- as.character(rangeZs$Chemical)
maxZs <- rangeZs[, .(zscore = max(zscore)), by = .(Idx, Chemical)]
maxZpos <- merge(maxZs, rangeZs, by = c("Idx", "Chemical", "zscore"))
peakPos <- maxZpos[, gp]
z.int.calc <- zscore.support.intervals(zdrop = 0.5) # was 0.5
zInts <- z.int.calc(supDTCpymlt)
findInts <- find.lod.dts(zInts, maxZpos)
#outerBounds <- findInts[, .(start = min(gp), end = max(gp)), by = Idx]
outerBounds <- findInts[, .(start = min(gp), end = max(gp), Chemical = as.character(unique(Chemical))), by = Idx]
flankIdxDT <- find.idx.range(supDTCpymlt, outerBounds)

comboZ <- flankIdxDT[, .(sumZs = sum(zscore)), by = .(Idx, gp)]
maxZsums <- comboZ[, .(sumZs = max(sumZs)), by = .(Idx)]
maxPos <- merge(maxZsums, comboZ, by = c("Idx", "sumZs"))
peakPosMax <- maxPos[, gp]
z.sum.int <- zscore.sums.support.intervals(zdrop = 0.75) # was 0.2
zSumsInt <- z.sum.int(comboZ)
findSumInts <- find.max.sum.zscores(zSumsInt, peakPosMax)
findSumIntsMrge <- merge(findSumInts, flankIdxDT, by = c("Idx", "gp"))[, c("zscore.x", "zscore2", "rightBound", "leftBound", "row") := NULL]
setnames(findSumIntsMrge, old = "zscore.y", new = "zscore")
findSumIntsMrge[, Chemical := as.character(Chemical)]
findSumIntsMrge[, Chemical := gsub("_", " ", Chemical)]

# ============================================================================
# Data wrangle the average frequency change of all haplotypes (for each chemical), with each chemical a different shape and the haplotypes colored consistently with previous plots.

chemShapes <- data.table(chems = unique(popDT$Chemical), shapes = c(0:6))
chemShapes2 <- chemShapes[, chems := gsub("18way_", "", chems)][, chems := gsub("_", " ", chems)]

avgDifsDT <- do.call(rbind, avgHapDifsDTs)
avgDifsDT[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
mergeAvgDifs <- merge(avgDifsDT, findSumIntsMrge, by = c("gp", "chr", "Chemical"))
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
        idxDT <- data.table(chr = allHaps$V1[1], pos = positions, gp = gpositions, Chemical = allHaps$V4[1], Idx = allHaps$V5[1], haps = hapVec2, freqDifs = allfreqs2, ID = idx$ID[1])
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
# Plotting genomewide zscores at positions where >= 2 founders are >=1.96 sd from the mean LOD score for each chemical

chroms <- as.character(as.roman(1:16))
start <- ch.bounds[seq(1, length(ch.bounds)-2, 2)]
end <- ch.bounds[seq(2, length(ch.bounds)-2, 2)]
chrDT <- data.table(start = start, end = end, chrBounds = rep("gray", 8))

gplot <- ggplot() + geom_point(data=supDTmlt, aes(gp, zscore, colour = Chemical), size = 2, alpha = 0.3) + ylab("zscore") + scale_colour_manual(values = colPal) + geom_hline(yintercept=1.96, linetype="dashed", color = "red") + geom_rect(data = chrDT, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.3) + theme_bw(base_size = 10) + theme(panel.grid = element_blank(), legend.position = "top") + scale_x_continuous(limits = c(0, ch.bounds[17]), breaks = c(g_l[1:16] + offsets[[2]][1:16]/2), labels = as.roman(1:16), expand = c(0, 0)) + xlab("")

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS18.pdf"), gplot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS18.png"), gplot, width = 8, height = 9, units = "in", dpi = 350)

# ============================================================================
# Plotting zscores for all genomic regions where may have pleiotropy (identified above) for each region separately (panel A), as well as the average haplotype change for each haplotype in each region separately (panel B)

idxAll <- do.call(rbind, idxLooperAll)
chemicals <- match(findSumIntsMrge$Chemical, chemShapes2$chems)
shapes <- chemShapes2$shapes[chemicals]
findSumIntsMrge[, Shape := shapes]
chrDT <- findSumIntsMrge[, .(chr = chr[1]), by = "ID"]

zPlot <-  ggplot(data=findSumIntsMrge, aes(`pos(kb)`, zscore, shape = Chemical, group = Chemical)) + geom_point(size = 2) + geom_line() + scale_shape_manual(values = setNames(findSumIntsMrge$Shape, findSumIntsMrge$Chemical)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), legend.position="none") + ylab("zscore") + xlab("position(kb)") + coord_cartesian(ylim = c(0, 4.5)) + facet_wrap(~ID, scales = "free_x", nrow = 2) + geom_text(data = chrDT, aes(x = -Inf, y = -Inf, hjust = -1.1, vjust = -1.25, label = paste0('chr', as.roman(chr))), inherit.aes = FALSE) + scale_x_continuous(labels = number_format(accuracy = 1))

shapeLeg <- data.table(shape = 0:6, chemical = c("cadmium chloride", "chlorpromazine", "diamide", "glacial acetic acid", "sodium chloride", "urea", "YPD")) 
colorLeg <- data.table(color.codes = unique(idxAll$color.codes), hapGrps = unique(idxAll$hapGrps))[order(hapGrps)]

avgHapPlot <- ggplot(data=idxAll, aes(`pos(kb)`, freqDifs, colour = hapGrps, shape = Chemical, group = interaction(finalGrp, Chemical))) + coord_cartesian(ylim = c(-1, 1)) + geom_point(size = 2) + geom_line() + scale_colour_manual(values=setNames(colorLeg$color.codes, colorLeg$hapGrps)) + scale_shape_manual(values = setNames(shapeLeg$shape, shapeLeg$chemical)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(0, "mm"), legend.spacing.y = unit(0, "mm"), plot.margin = unit(c(0, 0, 0, 0), "cm")) + ylab("mean haplotype change") + xlab("position (kb)") + facet_wrap(~ID, scales = "free_x", nrow = 2) + geom_text(data = chrDT, aes(x = -Inf, y = -Inf, hjust = -1.1, vjust = -1.25, label = paste0('chr', as.roman(chr))), inherit.aes = FALSE) + scale_x_continuous(labels = number_format(accuracy = 1))

legendPlot <- get_legend(avgHapPlot + guides(color = guide_legend(nrow = 2, override.aes = list(size = 2))) + guides(shape = guide_legend(nrow = 1, override.aes = list(size = 2))) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))
avgHapPlot <- avgHapPlot + theme(legend.position="none")

savePlot <- ggarrange(zPlot, avgHapPlot, ncol = 1, common.legend = TRUE, legend = "bottom", legend.grob = legendPlot, align = "v")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Pleiotropy_plots/FigureS19.pdf"), savePlot, width = 8, height = 10, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Pleiotropy_plots/FigureS19.png"), savePlot, width = 8, height = 10, units = "in", dpi = 350)

# ============================================================================
# Making table of genes under pleiotropic peaks

gff_file <- "saccharomyces_cerevisiae_R64-2-1_20150113.gff"

tags_to_keep <- c("X_element", "telomeric_repeat", "gene", "ARS", "long_terminal_repeat", "ncRNA_gene", "tRNA_gene", "snoRNA_gene", "centromere", "LTR_retrotransposon", "transposable_element_gene", "pseudogene", "Y_prime_element", "telomerase_RNA_gene", "snRNA_gene", "silent_mating_type_cassette_array", "mating_type_region", "intein_encoding_region", "rRNA_gene", "external_transcribed_spacer_region", "internal_transcribed_spacer_region", "non_transcribed_region", "origin_of_replication")

myGranges_plot <- gff_to_granges(gff_file, tags_to_keep )

peakTable <- pleio.peak.table(findSumIntsMrge, myGranges_plot)
fwrite(peakTable, paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/pleiotropy_peakTableRevised.txt"), sep = "\t", col.names = T)  

# ============================================================================
# Trouble-shooting