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
library(abind)
library(job)
library(ggdendro)
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

blah = function(v){
    F = unlist(strsplit(as.character(v[1]),";"))
    fD = unlist(strsplit(as.character(v[2]),";"))
    oo = order(as.numeric(fD),decreasing=TRUE)
    oF = F[oo]
    ofD = fD[oo]
    # most changed, 2nd most, 3rd most, 10th most
    c(oF[1],ofD[1],oF[2],ofD[2],oF[3],ofD[3],oF[10],ofD[10])
}

# ============================================================================
# Load data needed for downstream analyses

founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)
indHaps <- read.table("Individual_haplotype_differences_55_selected_populations.txt",header=TRUE)
indHapsDT <- as.data.table(indHaps)
avgHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Averaged_and_sd_tx_tables/", analysisType = "hap_freqs_avg_sd_difs", samplePattern = ".*")
avgHapDifsDTsCpy <- copy(avgHapDifsDTs)


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
gray <- col2rgb('grey50')
chromBarCol <- rgb(gray[1],gray[2],gray[3], maxColorValue = 255, alpha = 100)
colorCodes = c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
hx = sapply(strsplit(colorCodes, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
hx2 = hx[-c(5,21,24)]   
treatmentAbbrv <- fread("chemAbbrev.txt")
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

names(colPal) <- plottingFactors
treatmentAbbrv <- fread("chemAbbrev.txt")
names(colPal) <- treatmentAbbrv$abbr
colorPal <- colPal[order(names(colPal))]
plottingFactors <- read.table("all_chems_list.txt", header = T, sep = "\t")
plottingFactors2 <- plottingFactors$Chemical
plottingFactors3 <- plottingFactors2[c(2,1,4,3,5,6,7,8,9,12,10,11,13,14,15,16)]
names(colorPal) <- plottingFactors3

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
topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))

# ============================================================================
# Find the two most increased haplotypes per position.

idxLooper <- lapply(avgHapDifsDTs, function(chem)  {
    tic()
    freqDifs <- chem[, tstrsplit(avgDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
    baseFreqs <- chem[, tstrsplit(baseFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
    avgFreqs <- as.data.frame(chem[, tstrsplit(avgFreq, split = ";", type.convert = TRUE, fixed = TRUE)])
    collHaps <- as.data.frame(chem[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
    baseFreqsDF <- as.data.frame(baseFreqs)
    maxVals <- apply(freqDifs, 1, max, na.rm = T) ## only looking at increasing values
    maxIdx <- apply(freqDifs, 1, which.max)
    maxHap <- collHaps[cbind(seq_along(maxIdx), maxIdx)]
    baseHap <- baseFreqsDF[cbind(seq_along(maxIdx), maxIdx)]
    evFreq <- avgFreqs[cbind(seq_along(maxIdx), maxIdx)]
    chem[, c("maxHap", "maxChange", "baseHap", "evFreq") := .(maxHap, maxVals, baseHap, evFreq)]
    nxtMax <- apply(freqDifs, 1, function(x) rev(sort(x))[2])
    nxtMaxIdx <- apply(freqDifs, 1, function(x) which(x == rev(sort(x))[2])[1])
    maxHap2 <- collHaps[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    baseHap2 <- baseFreqsDF[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    evFreq2 <- avgFreqs[cbind(seq_along(nxtMaxIdx), nxtMaxIdx)]
    chem[, c("maxHap2", "maxChange2", "baseHap2", "evFreq2") := .(maxHap2, nxtMax, baseHap2, evFreq2)]
    idx <- chem[, c("chr", "pos", "gp", "Chemical", "maxHap", "maxChange", "maxHap2", "maxChange2", "baseHap", "baseHap2", "evFreq", "evFreq2")]
    toc()
    #idx[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
    idx
} )

rep1Freqs <- lapply(idxLooper, function(diffs) {
    tic()
    indDT <- indHapsDT[Chemical == diffs$Chemical[1]]
    indDT <- indDT[Replicate == unique(Replicate)[1]]
    hapFreqs <- as.data.frame(indDT[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)])
    collHaps <- as.data.frame(indDT[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
    findMax <- unlist(lapply(1:nrow(collHaps), function(hap) {
        which(collHaps[hap,] == diffs$maxHap[hap])
    }))
    maxFreqs <- hapFreqs[cbind(seq_along(findMax), findMax)]
    diffs$rep1freqs <- maxFreqs
    indDT2 <- indHapsDT[Chemical == diffs$Chemical[1] & Replicate != unique(Replicate)[1]]
    hapFreqs2 <- as.data.frame(indDT2[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)])
    findAllMax <- rep(findMax, indDT2[,.N]/indDT[,.N])
    maxFreqs2 <- hapFreqs2[cbind(seq_along(findAllMax), findAllMax)]
    allGP <- rep(diffs$gp, indDT2[,.N]/indDT[,.N])
    allDT <- data.table(gp = allGP, maxFreqs = maxFreqs2)
    allDTavg <- allDT[, .(avgFreq = mean(maxFreqs)), by = gp]
    diffs$restOfReps <- allDTavg$avgFreq
    diffs[, absDiffs := abs(rep1freqs - restOfReps)]
    toc()
    diffs
})
chemMrge <- do.call(rbind, rep1Freqs)
fwrite(chemMrge, paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/MIH_deviation_DT.txt"), sep = "\t", col.names = T)  

# ============================================================================
# Plot the absolute difference of the MIH b/w the first replicate and the average of the remaining replicates for each treatment genomewide

mihDevDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/MIH_deviation_DT.txt"), sep = "\t", header = T)

filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Ind_reps_hap_dev_plots/")
gray <- col2rgb('grey50')
chromBarCol <- rgb(gray[1],gray[2],gray[3], maxColorValue = 255, alpha = 100)
txDTs <- split(mihDevDT, mihDevDT$Chemical)

job::job({
textSize <- 1.5
ylabel <- "Absolute haplotype deviation of the MIH"
plotCounter <- 0
pdf(paste0(filePath, "SEE01_all_chems_MIH_hap_devs.pdf"), width = 8, height = 10)
#png(paste0(filePath, "SEE01_all_chems_MIH_hap_devs.png"), width = 8, height = 10, units = 'in', res = 350)
layout(matrix(1:length(txDTs)))
par(mar = c(0.7,2.5,0,6), mgp = c(1.5,0.7,0), oma = c(2,2.5,2,2))
plotting <- lapply(txDTs, function(DT) {
    plotCounter <<- plotCounter + 1
    tic()
    chemUsing <- gsub("18way_|_12", "", DT$Chemical[1])
    chemUsing <- gsub("_", " ", chemUsing)
    print(chemUsing)
    flush.console()
    ymax <- 1
    ymin <- 0
    DT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(gp),max(gp)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
    DT[, box()]
    #DT[, polygon(x = c(gp, rev(gp)), y = c(lowerCI, rev(upperCI)), col =  adjustcolor("dodgerblue", alpha.f = 0.25), border = NA)]
    DT[, points(gp, absDiffs, col = "black", bty = "n", pch = 16, cex = 0.5)]
    chromSplit <- split(DT, DT$chr)
    smoothing <- lapply(chromSplit, function(ch) {
        ch[, lines(ksmooth(gp, absDiffs, kernel = "normal", bandwidth = 100000), col = "red", lwd = 1)]
    } )
    lapply(seq(1, length(ch.bounds)-2, 2), function(y) {
        rect(ch.bounds[y], ymin, ch.bounds[y+1], ymax, col = chromBarCol, border = NA, xpd = TRUE) 
    } )
    DT[, text(x = 10500000, y = 0.9, labels = chemUsing, cex = textSize)]
    DT[, axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1), las = 1, cex.axis = 0.8)]
    #DT[, mtext(text = yTickLab, side = 2, line = 0.5, cex = 0.6, las = 1)]
    #if(plotCounter == 1) {
        #DT[, title(main = mainTitle, cex.main = 1, line = 0.5, xpd = NA)]}
    if(length(txDTs) %% 2 == 1 & plotCounter == length(txDTs)/2 + 0.5) {
        DT[, mtext(text = ylabel, side = 2, line = 3.25, cex = textSize)]} else if(length(txDTs) %% 2 == 0 & plotCounter == round(length(txDTs)/2, 0)) {DT[, text(grconvertX(-0.07, "npc", "user"), grconvertY(0.02, "npc", "user"), labels = ylabel, xpd = NA, cex = textSize, srt = 90)]}
    if(plotCounter == length(txDTs)) {DT[, axis(1, at=g_l[1:16] + offsets[[2]][1:16]/2, labels= as.roman(1:16), las = 1, cex.axis = 1)]}
    toc()
} )
dev.off()
}, import = c(txDTs, filePath, ch.bounds, chromBarCol, g_l, offsets), packages = c("ggplot2", "data.table", "GGally", "ggpubr", "tictoc"))

# ============================================================================
# Make a table of the frequency change of the MIH for all replicates at chromosome 2 for a each chemical to zoom-in on how replicable adaptation is 
#chemShapes <- data.table(chems = unique(popDT$Chemical), shapes = c(0:6))
#chemShapes2 <- chemShapes[, chems := gsub("18way_", "", chems)][, chems := gsub("_", " ", chems)]
#chromosome <- 5

repFreqs <- lapply(idxLooper, function(diffs) {
    tic()
    indDT <- indHapsDT[Chemical == diffs$Chemical[1]] 
    indDT[, gp := rep(diffs$gp, uniqueN(Replicate))]
    #indDT <- indDT[chr == chromosome]
    #diffs <- diffs[chr == chromosome]
    hapFreqs <- as.data.frame(indDT[, tstrsplit(freqDifs, split = ";", type.convert = TRUE, fixed = TRUE)])
    collHaps <- as.data.frame(indDT[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
    diffreps <- rep(diffs$maxHap, indDT[,uniqueN(Replicate)])
    findMax <- unlist(lapply(1:nrow(collHaps), function(hap) {
        which(collHaps[hap,] == diffreps[hap])
    }))
    maxFreqs <- hapFreqs[cbind(seq_along(findMax), findMax)]
    indDT$maxFreqChange <- maxFreqs
    indDT$maxHap <- diffreps
    hapColors <- match(indDT$maxHap, topHaps)
    hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
    Colors <- mycols[hapColors]
    indDT[, "color.codes":= Colors][, "pos(kb)" := pos/1000]
    chem <- indDT[, repGrps := Replicate][order(Replicate, gp)][, c("difCol") := .(c(1000, diff(gp))), by = "Replicate"][, row := .I]
    chem[, grpCol := paste0(Replicate, "_", difCol)]
    chem2 <- chem[, grp := cumsum(c(TRUE, diff(gp)!=1000)), by = Chemical][order(row)] ### troubleshoot have two different haps in same group (gp 4207890)
    chem2[, grp := cumsum(c(TRUE, diff(gp)!=1000)), by = Replicate]
    chem2[, grp2 := rleid(Replicate)]
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
idxDTs <- do.call(rbind, repFreqs)
idxGrps <- rle(idxDTs$grp)
idxGrps$values <- 1:length(idxGrps$values)
idxDTs$grp <- rep(idxGrps$values, idxGrps$lengths)

fwrite(idxDTs, paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/MIH_allReps_DT.txt"), sep = "\t", col.names = T)  

# ============================================================================
# Plot the MIH for all replicates at chromosome 2 for a each chemical to zoom-in on how replicable adaptation is 

mihAllDT <- idxDTs
#mihAllDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/MIH_allReps_DT.txt"), sep = "\t", header = T)
mihAllDT[, Chemical := gsub("18way_|_12", "", Chemical)]
mihAllDT[, Chemical := gsub("_", " ", Chemical)]
colorLeg <- data.table(color.codes = mycols, haps = c(topHaps, "other"))
#colorLeg <- data.table(color.codes = unique(mihAllDT$color.codes), repGrps = unique(mihAllDT$repGrps))[order(repGrps)]

avgHapPlot <- ggplot(data=mihAllDT[Chemical == "cadmium chloride"], aes(`pos(kb)`, maxFreqChange, colour = maxHap, group = interaction(finalGrp, Replicate))) + coord_cartesian(ylim = c(-.5, 1)) + geom_point(size = 0.5) + geom_line(color="lightgrey", size=0.25) + scale_colour_manual(values=setNames(colorLeg$color.codes, colorLeg$haps)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(0, "mm"), legend.spacing.y = unit(0, "mm"), plot.margin = unit(c(0, 0, 0, 0), "cm")) + ylab("Haplotype frequency change") + xlab("position on chrV (kb)") + facet_wrap(~Chemical, scales = "free_x", nrow = 7) + scale_x_continuous(labels = number_format(accuracy = 1)) + theme(legend.position = "right") + guides(color = guide_legend(ncol = 2, override.aes = list(size = 1)))

#legendPlot <- get_legend(avgHapPlot + guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))
#avgHapPlot <- avgHapPlot + theme(legend.position="none")
#savePlot <- ggarrange(avgHapPlot, ncol = 1, common.legend = TRUE, legend = "right", legend.grob = legendPlot, align = "v")

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Pleiotropy_plots/FigureSx_replicable_adaptation.pdf"), avgHapPlot, width = 8, height = 4, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Pleiotropy_plots/FigureSx_replicable_adaptation.png"), avgHapPlot, width = 8, height = 4, units = "in", dpi = 350)

# ============================================================================
# Plot the correlation in average absolute difference between replicate 1 and the remaining replicates for the MIH with the average heterozygosity 
mihDevDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/MIH_deviation_DT.txt"), sep = "\t", header = T)

hapDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"), header = T)
hapDT$treatment <- gsub("^18way_|_([0-1]?[0-9])$", "", hapDT$chemical)
avgHets <- hapDT[, mean(het), by = chemical]
setnames(avgHets, c("chemWeek", "avgHet"))
avgHets$id <- unlist(lapply(avgHets$chemWeek, function(x) {
    splitting <- strsplit(x,"_")[[1]]
    if(length(splitting) == 3) {return(splitting[2])} else{
        return(paste(splitting[2:3], collapse = "_"))}
} ) )

avgChemCorDFs <- lapply(unique(chemNames), function(x) {
    chemDF <- corJob[grepl(x, rownames(corJob)), grepl(x, names(corJob))]
    avgCor <- mean(chemDF[upper.tri(chemDF)])
    avgCorDF <- data.frame(chemical = x, avgCor = avgCor, stringsAsFactors = FALSE)
    avgCorDF
} )
avgDevDT <- mihDevDT[, .(avgDiff = mean(absDiffs)), by = "Chemical"][, chemical := gsub("18way_|_12", "", Chemical)]

avgHets[, chemical := gsub("18way_|_12", "", chemWeek)][, chemWeek := NULL]
avgCorHetDF <- merge(avgHets, avgDevDT , by = "chemical")
avgCorHetDF[chemical %like% "YPD", c("chemical", "id") := .('ypd', 'ypd')]
avgCorHetDF <- avgCorHetDF[order(avgCorHetDF$chemical),]
avgCorHetDF <- na.omit(avgCorHetDF)
avgCorHetDF$chemical[7] <- "ypd"

#chemsUsing <- gsub("_", " ", avgCorHetDF$chemical)
#chemsUsing[chemsUsing == "YPD"] <- "ypd"
avgCorHetDF$chemical <- gsub("_", " ", avgCorHetDF$chemical)


gplot <- ggplot(data=avgCorHetDF, aes(avgHet, avgDiff)) + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(aes(label = ..rr.label..)) + geom_point(aes(colour = chemical), size = 2) + xlab("average within-treatment heterozygosity") + ylab("average within-treatment deviation") + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + scale_colour_manual(values = colorPal[names(colorPal) %in% avgCorHetDF$chemical])

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS13_revised.pdf"), gplot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS13_revised.png"), gplot, width = 8, height = 9, units = "in", dpi = 350)

# ============================================================================
# Find the average number of haplotypes detected across all 55 selected populations per position


collHaps <- as.data.frame(indHapsDT[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)])
names(collHaps) <- paste0("V", 2:18)
indHapsDT[, chrPos := paste0(chr, "_", pos)]
collHaps2 <- cbind(indHapsDT$chrPos, collHaps)
names(collHaps2)[1] <- "chrPos"
detectHaps <- unlist(lapply(1:nrow(collHaps2), function(x) {
    length(which(collHaps2[x,2:18] > 0.0003))
}))

reps <- length(unique(indHapsDT$id))
df2 <- as.data.frame(matrix(detectHaps, byrow=TRUE, ncol = 55))
avgHapsDetected <- rowMeans(df2)
avgHapsPerSite <- mean(avgHapsDetected)

fwrite(as.data.table(df2), "haps_detected_per_position_55_replicates.txt", col.names = T, sep = "\t")

indHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_sum_of_squared_diffs_tables/", analysisType = "haps_sq_diffs", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]_")
dt1 <- indHapDifsDTs[[1]]
collHapFreqs <- as.data.frame(dt1[, tstrsplit(baseRawFreqsCol, split = ";", type.convert = TRUE, fixed = TRUE)])
names(collHapFreqs) <- paste0("V", 2:18)
collHapFreqs2 <- cbind(dt1$gp, collHapFreqs)
names(collHapFreqs2)[1] <- "gp"
detectedHaps <- unlist(lapply(1:nrow(collHapFreqs2), function(x) {
    length(which(collHapFreqs2[x,2:18] > 0.005))
}))
avgHapsPerSite <- mean(detectedHaps)




