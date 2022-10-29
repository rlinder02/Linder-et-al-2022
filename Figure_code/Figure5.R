library(tictoc)
library(data.table)
library(qqman)
library(rtracklayer)
library(RColorBrewer)
library(Gviz)
library(scales)
library(R.utils)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggbeeswarm)
library(viridis)
#require(svMisc)
library(GenomicRanges)

source('formatting/Haplotype_file_splitter.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')
source('calculating/Calculating_test_statistics.R')

# ============================================================================
# Set global options

defDir <- getwd()
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))
options("scipen"=100, "digits"=10)
# ============================================================================
# Custom functions

pairwise.founder.compare <- function(founder1, founder2, maxDiff) {
    function(DT) {
        newCol <- paste0(founder1, "-", founder2)
        DT[, newCol := get(paste0(founder1, ".difs")) - get(paste0(founder2, ".difs"))]
        DT[newCol > -maxDiff & newCol < maxDiff]
    }
}

findSigRows <- function(vec, condition, value) {
    function(x) {
        findSigs <- sort(unique(unlist(lapply(vec, function(hap) {
            which(condition(x[, get(hap)], value))
        }))))
        findSigs
    }        
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

plottingWindows <- function(peakPos, DT, flanking) {
    counter <- 0
    findFlanks <- lapply(peakPos, function(peak) {
        counter <<- counter + 1
        DT[(gp >= peak - flanking) & (gp <= peak + flanking)][, Idx := counter]
    } )
    flanksDT <- do.call(rbind, findFlanks)
    flanksDT
} 

add.alpha <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
              rgb(x[1], x[2], x[3], alpha=alpha))  
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

rep.list <- function(object, repObject) {
    rep(list(object), length(repObject))
}

peak.table <- function(lodWin, myGranges2, numPeaks) {
    tic('Total')
    findTopPeaks <- lodWin[significant == 1][order(-LOD, )][1:numPeaks][, Idx]
    topPeaksDT <- lodWin[Idx %in% findTopPeaks]
    splitIdxs <- split(topPeaksDT, topPeaksDT$Idx)
    genomeCoordinates <- lapply(splitIdxs, function(win) {
        chrName <- as.roman(win$chr[1])
        chr <- paste0('chr', chrName)
        start <- min(win$pos)
        end <- max(win$pos)
        interval <- end-start
        chem <- gsub("18way_", "", win$Chemical[1])
        coordsDF <- data.frame(chrom = chr, start = start, end = end)
        coordsGR <- makeGRangesFromDataFrame(coordsDF)
        overlappingRegion <- subsetByOverlaps(myGranges2, coordsGR)
        genes <- overlappingRegion$symbol
        geneDT <- data.table(chemical = chem, chr = chr, start_pos = start, end_pos = end, interval_length = interval/1000, genes = I(list(genes)))
        setnames(geneDT, "interval_length", "interval(kb)")
    } )
    table1 <- do.call(rbind, genomeCoordinates)
    toc()
    table1		
}

plotting_windows_gviz <- function(plotAllWins, avgDifs, folder, statistic, snpFreqs) {
    #topHaps <- topHapsDT[chemWeek == plotAllWins$chemWeek[1], topHapCombos]
    #topHaps <- strsplit(topHaps, ";")[[1]]
    snpFreqs2 <- dcast(snpFreqs, chr + POS + gp + id ~ founder, value.var = "control.difs")
    snpFreqs2 <- as.data.frame(snpFreqs2, stringsAsFactors = FALSE)
    counter <- 0
    plotWinsList <- split(plotAllWins, plotAllWins$Idx)
    plotting <- lapply(plotWinsList, function(x) {
        counter <<- counter + 1
        print(counter)
        flush.console()
        treatment <- x$id[1]
        if(x$chr[1] == 'chrmt') {chr_name <- 'M'} else {
            chr_name <- as.roman(x$chr[1])}
        fname <- paste0(treatment, "_", "chr", chr_name, "_", min(x$pos), "..", max(x$pos))
        print(fname)
        flush.console()
        chr <- chr_name
        start <- min(x$pos)
        end <- max(x$pos)
        interval <- nrow(x)
        sigHapDifs <- avgDifs[chr == x$chr[1] & pos >= start & pos <= end]
        sigDifs <-sigHapDifs[, tstrsplit(avgDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
        sigHaps <-sigHapDifs[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
        topSix <- lapply(1:nrow(sigDifs), function(x) {
            topVals <- sigDifs[x, order(abs(as.numeric(sigDifs[x])), decreasing = TRUE)[1:6]]
            topFounders <- sigHaps[x, ..topVals]
            topValues <- sigDifs[x, ..topVals]
            foundVals <- cbind(sigHapDifs[x, pos], topFounders, topValues)
            setnames(foundVals, c("pos", paste0("V", 1:(ncol(foundVals)-1))))
        } )
        topBind <- do.call(rbind, topSix)
        xstart <- sort(rep(topBind$pos - 400, 6))
        xend <- sort(rep(topBind$pos + 400, 6))
        ystart <- as.numeric(unlist(t(topBind[, V7:V12])))
        yend <- ystart
        hapVec <- as.character(unlist(t(topBind[, V1:V6])))
        hapColors <- match(hapVec, topHaps)
        hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
        hapThick <- lengths(regmatches(hapVec, gregexpr("[A-B][0-9]", hapVec))) + 2
        subsettingSnps <- snpFreqs2[snpFreqs2$gp >= min(x$gp) & snpFreqs2$gp <= max(x$gp),]
        snpVals <- subsettingSnps[, 5:ncol(subsettingSnps)]
        reorderSnpVals <- unlist(lapply(founderNames, function(reord) {
            grep(reord, names(snpVals))
        }))
        snpVals2 <- snpVals[c(reorderSnpVals, 1)]
        if(nrow(subsettingSnps) == 0) {snpY <- c(0,0); grSnps <- GRanges(seqnames = paste0('chr', chr), ranges = IRanges(start = x$pos, end = x$pos))} else {
            maxSnpY <- max(subsettingSnps[,5:ncol(subsettingSnps)], na.rm = TRUE); snpY <- pretty(c(0,maxSnpY)); grSnps <- GRanges(seqnames = paste0('chr', chr), ranges = IRanges(start = subsettingSnps$POS, end = subsettingSnps$POS)); mcols(grSnps) <- snpVals2}
        snpYLims <- c(snpY[1], snpY[length(snpY)])
        mycolorsSnps <- match(names(snpVals2), topHaps)
        mycolorsSnps[is.na(mycolorsSnps)] <- (length(topHaps) + 1)
        maxHapY <- 1
        minHapY <- -1
        hapY <- c(-1, -0.5, 0, 0.5, 1)
        hapYLims <- c(-1, 1)
        maxStatY <- x[, max(get(statistic))]
        minStatY <- x[, min(get(statistic))]
        statY <- pretty(c(minStatY, maxStatY))
        statYLims <- c(statY[1],statY[length(statY)])
        maxSumSq <- max(sigHapDifs$sumSqDifs)
        minSumSq <- 0
        sumSqY <- pretty(c(minSumSq, maxSumSq))
        sumSqYLims <- c(sumSqY[1],sumSqY[length(sumSqY)])
        emptyDF <- as.data.frame(matrix(data = -2, ncol = length(topHaps)+1, nrow = nrow(x)))
        names(emptyDF) <- c(topHaps, "other")
        emptyDF[] <- lapply(emptyDF, as.numeric)
        grHaps <- GRanges(seqnames = paste0('chr', chr), ranges = IRanges(start = x$pos, end = x$pos))
        mcols(grHaps) <- emptyDF
        groupingHaps <- names(mcols(grHaps))
        groupingHaps <- factor(groupingHaps, levels = groupingHaps)
        groupingSnps <- names(mcols(grSnps))
        groupingSnps <- factor(groupingSnps, levels = groupingSnps)
        grStat <- GRanges(seqnames = paste0('chr', chr), ranges = IRanges(start = x$pos, end = x$pos))
        mcols(grStat) <- x[, get(statistic)]
        grSumSq <- GRanges(seqnames = paste0('chr', chr), ranges = IRanges(start = x$pos, end = x$pos))
        mcols(grSumSq) <- sigHapDifs[, "sumSqDifs"]
        dtrackSnps <- DataTrack(range = grSnps, genome = "sacCer3", groups = groupingSnps, col = mycols[mycolorsSnps], name = "Absolute SNP \nDifferences", ylim = snpYLims, legend = FALSE)
        hapYlims <- c(min(ystart) - 0.25*abs(min(ystart)), max(ystart) + 0.25*abs(max(ystart))) 
        dtrackHaps <- DataTrack(range = grHaps, genome = "sacCer3", name = "Average haplotype \nchange", ylim = hapYlims)
        cTrack <- CustomTrack(plottingFunction = function(GdObject, prepare) {
            if(!prepare) pushViewport(dataViewport(xData = c(xstart, xend), yData = ystart, xscale = NULL, yscale = hapYlims))
            if(!prepare) grid.segments(x0 = xstart, y0 = ystart, x1 = xend, y1 = yend, default.units = "native", arrow = NULL, name = NULL, gp = gpar(col = mycols[hapColors], lty = "solid", lwd = hapThick), draw = TRUE)
            #if(!prepare) grid.legend(c(topHaps, "other"), nrow = 4, ncol = 6, pch = 16, hgap = 0.5, vgap = 0.5, gp = gpar(col = mycols, cex = 0.5), draw = TRUE)
            return(invisible(GdObject))
        } )
        ctrackLegend <- CustomTrack(plottingFunction = function(GdObject, prepare) {
            if(!prepare) grid.legend(c(topHaps, "other"), nrow = 4, ncol = 6, pch = 16, hgap = 1, vgap = 0.5, gp = gpar(col = mycols, cex = 0.6), draw = TRUE)
            return(invisible(GdObject))
        } )
        displayPars(ctrackLegend) <- list(showTitle = FALSE, showAxis = FALSE)
        ot <- OverlayTrack(trackList=list(dtrackHaps, cTrack))
        dtrackStat <- DataTrack(range = grStat, genome = "sacCer3", name = "-log10(p)", type = "b", col = "black", frame = FALSE, lwd = 3, jitter.x = FALSE, cex = 1, ylim = statYLims, lty = 1)
        dtrackSumSq <- DataTrack(range = grSumSq, genome = "sacCer3", name = "Sum of \nsquared \nhap change", type = "b", col = "red", frame = FALSE, lwd = 3, jitter.x = FALSE, cex = 1, ylim = sumSqYLims, lty = 1)
        grtrack <- GeneRegionTrack(myGranges_plot, chromosome=paste0('chr',chr), start=start, end=end, showId=TRUE, shape='fixedArrow', thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5", "miRNA", "lincRNA"), showFeatureId=TRUE, showTitle = TRUE, name = 'Genes', just.group = "below", geneSymbol = TRUE, feature = myGranges_plot$feature, featureAnnotation = unique(myGranges_plot$feature), arrowHeadWidth=10)
        if(as.character(chr) == 'mt') {itrack_chr <- 'chrM'} else {itrack_chr <- paste0('chr', chr)}
        gtrack <- GenomeAxisTrack(showTitle = TRUE, name = itrack_chr, labelPos = "alternating", littleTicks = FALSE, rotation.title = 0)
        if(nrow(subsettingSnps) == 0) {
            items_to_plot <- list(dtrackStat, dtrackSumSq, gtrack, grtrack, ctrackLegend, ot)
            relSizes <- c(3,3,1,3,2,4)} else{
                items_to_plot <- list(dtrackStat, dtrackSumSq, gtrack, grtrack, dtrackSnps, ctrackLegend, ot)
                relSizes <- c(3,3,1,3,3,2,4)
            }
        pdf(paste0(folder, fname,".pdf"), height = 12, width = 8)
        plotTracks(items_to_plot, from = start, to = end, cex.axis = 1, cex.title = 1, col.title = "black", col.axis = "black", chromosome = paste0('chr',chr), title.width = 1, sizes = relSizes)
        dev.off()
    } )
}

# ============================================================================
# Load data

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
avgHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Averaged_and_sd_tx_tables/", analysisType = "hap_freqs_avg_sd_difs", samplePattern = ".*")
indHapDifsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_sum_of_squared_diffs_tables/", analysisType = "haps_sq_diffs", samplePattern = "^SEE12B02")
indHapDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_het_diffs_tables/", analysisType = "haps_het_diffs", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]_")
indHapDTsUsing <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/All_reps_tx_tables/", analysisType = "hap_freqs", samplePattern = "^18way")
allBind <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_All_Reps_Hets.txt"))
topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))
hetDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"), header = T)

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
gray <- col2rgb('grey50')
chromBarCol <- rgb(gray[1],gray[2],gray[3], maxColorValue = 255, alpha = 100)
colorCodes = c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
hx = sapply(strsplit(colorCodes, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
hx2 = hx[-c(5,21,24)]   
mycols = c(hx2,"#000000")

# ============================================================================
# Find significant peaks from haplotype-corrected LOD scores

colNames <- rep(list('LOD'), length(indLODDTs))
#threshold <- rep(list(25), length(indLODDTs))
threshold <- rep(list(50), length(indLODDTs))
allPeaks <- Map(inflect, indLODDTs, colNames, threshold)
peakCol <- Map(add.sig.column, indLODDTs, allPeaks, rep.list(5, indLODDTs))

# ============================================================================
# Plotting haplotype-adjusted LOD scores for the averaged treatments vs base genomewide

chroms <- as.character(as.roman(1:16))
filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Avg_founder_hap_adjusted_smoothed_LOD/")
counter <- 0
gwidePlot <- lapply(indLODDTs, function(DT) {
    counter <<- counter + 1
    fname <- DT$chemWeek[1]
    print(fname)
    flush.console()
    ymax <- max(DT$LOD) + 0.02*max(DT$LOD)
    ymin <- 0
    yTicks <- pretty(ymin:ymax)
    pdf(paste0(filePath, fname, "_vs_collapsed_adjusted_LOD_conservative.pdf"),width = 6, height = 2)
    par(mar = c(2,3,0,0), mgp = c(1,0.5,0))
    DT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(gp),max(gp)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
    DT[, lines(gp, LOD, col = "black", bty = "n", lwd = 1)]
    abline(h = 5, col = "red", lty = 3)
    DT[, box()]
    lapply(seq(1, length(ch.bounds), 2), function(y) {
        rect(ch.bounds[y], ymin, ch.bounds[y+1], ymax, col = chromBarCol, border = NA, xpd = TRUE) 
    } )
    DT[, axis(1, at=g_l[1:17] + offsets[[2]][1:17]/2, labels= as.roman(c(1:16, 'm')), las = 1, cex.axis = 0.6)]
    DT[, axis(2, at=yTicks, labels=yTicks, las = 1, cex.axis = 0.6)]
    DT[, mtext(text = "-log10(p)", side = 2, line = 2, cex = 0.75, las = 0)]
    peaks <- peakCol[[counter]]
    peaks <- peaks[significant == 1]
    points(peaks$gp, peaks$LOD, pch = 16, col = "red", cex = 0.5)
    graphics.off()
} )

# ============================================================================
# Plotting haplotype-adjusted LOD scores for the averaged treatments vs base for cadmium chloride chrXVI

filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Avg_founder_hap_adjusted_smoothed_LOD/")
allLODDTs <- do.call(rbind, indLODDTs)
txDTs <- split(allLODDTs, allLODDTs$Chemical)
chrom <- 16
DT <- txDTs[[1]]
DT <- DT[chr == chrom]
peaks <- peakCol[[1]]
peaks <- peaks[chr == 16]
peaks <- peaks[significant == 1]

fname <- DT$Chemical[1]
print(fname)
flush.console()
ymax <- max(DT$LOD) + 0.02*max(DT$LOD)
ymin <- 0
xmin <- DT[, min(pos)]
xmax <- DT[, max(pos)]
xticks <- pretty(xmin:xmax)
yTicks <- pretty(ymin:ymax)
pdf(paste0(filePath, fname, "_chr", as.roman(chrom),  "_hap_adjusted_LOD.pdf"),width = 6, height = 2)
par(mar = c(2,3,0,0), mgp = c(1,0.5,0))
DT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(xmin,xmax), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
DT[, lines(pos, LOD, col = "black", bty = "n", lwd = 2)]
points(peaks$pos, peaks$LOD, pch = 16, col = "red", cex = 0.5)
abline(h = 5, col = "red", lty = 3)
DT[, box()]
DT[, axis(1, at=xticks, labels= xticks/1000, las = 1, cex.axis = 0.6)]
DT[, axis(2, at=yTicks, labels=yTicks, las = 1, cex.axis = 0.6)]
DT[, mtext(text = "-log10(p)", side = 2, line = 2, cex = 0.75, las = 0)]
DT[, text(x = 890000, y = 65, labels = paste0("chr", as.roman(chrom)), cex = 0.75)]
DT[, mtext(text = "position(kb)", side = 1, line = 1, cex = 0.75, las = 0)]
graphics.off()

# ============================================================================
# Manually specify peaks from haplotype-corrected LOD scores

allLODDTs <- do.call(rbind, indLODDTs)
cadRegion1 <- allLODDTs[chemWeek == "18way_cadmium_chloride_12"][chr == "16"][pos >= 160000 & pos <= 175000][, c("significant", "Idx") := .(0, 1)]
#diaRegion1 <- allLODDTs[chemWeek == "18way_diamide_12"][chr == "15"][pos >= 260000 & pos <= 270000][, c("significant", "Idx") := .(0, 2)]
#naclRegion1 <- allLODDTs[chemWeek == "18way_sodium_chloride_12"][chr == "4"][pos >= 530000 & pos <= 535000][, c("significant", "Idx") := .(0, 3)]
#chlorpRegion1 <- allLODDTs[chemWeek == "18way_chlorpromazine_12"][chr == "1"][pos >= 35000 & pos <= 45000][, c("significant", "Idx") := .(0, 4)]
#chlorpRegion2 <- allLODDTs[chemWeek == "18way_chlorpromazine_12"][chr == "2"][pos >= 330000 & pos <= 340000][, c("significant", "Idx") := .(0, 5)]
#chlorpRegion3 <- allLODDTs[chemWeek == "18way_chlorpromazine_12"][chr == "12"][pos >= 805000 & pos <= 820000][, c("significant", "Idx") := .(0, 6)]
#sigRegions <- rbind(cadRegion1, diaRegion1, naclRegion1)
sigRegions <- cadRegion1
# ============================================================================
# Custom gviz plots (user-defined intervals from sigRegions defined above)

topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))

gff_file <- "saccharomyces_cerevisiae_R64-2-1_20150113.gff"

tags_to_keep <- c("X_element", "telomeric_repeat", "gene", "ARS", "long_terminal_repeat", "ncRNA_gene", "tRNA_gene", "snoRNA_gene", "centromere", "LTR_retrotransposon", "transposable_element_gene", "pseudogene", "Y_prime_element", "telomerase_RNA_gene", "snRNA_gene", "silent_mating_type_cassette_array", "mating_type_region", "intein_encoding_region", "rRNA_gene", "external_transcribed_spacer_region", "internal_transcribed_spacer_region", "non_transcribed_region", "origin_of_replication")

myGranges_plot <- gff_to_granges(gff_file, tags_to_keep )

matchIds <- unlist(lapply(unique(sigRegions$chemWeek), function(x) {
    counter <- 0
    findIdx <- unlist(lapply(avgHapDifsDTs, function(y) {
        counter <<- counter + 1
        if(y$chemWeek[1] == x) {return(counter)}
    } ) )
    return(findIdx)
} ) )

matchSnpIds <- unlist(lapply(unique(sigRegions$id), function(x) {
    counter <- 0
    findIdx <- unlist(lapply(avgSnpDifsDTs2, function(y) {
        counter <<- counter + 1
        ident <- strsplit(y$id[1], "-")[[1]][1]
        if(grepl(ident, x)) {return(counter)}
    } ) )
    return(findIdx)
} ) )

avgHapDifsDTs2 <- avgHapDifsDTs[matchIds]
avgSnpDifsDTs3 <- avgSnpDifsDTs2[matchSnpIds]

plotAllWins <- split(sigRegions, sigRegions$chemWeek)
avgDifs <- avgHapDifsDTs2
statistic <- rep.list("LOD", plotAllWins)
folder <- rep.list(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/LOD_avg_squared_difs_vs_base_collapsed_plot_windows/"), plotAllWins)
snpFreqs <- avgSnpDifsDTs3

plotLODHapSnpDifs <- Map(plotting_windows_gviz, plotAllWins, avgDifs, folder, statistic, snpFreqs)

# ============================================================================
# Trouble-shooting
