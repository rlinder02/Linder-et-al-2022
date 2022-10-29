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
library(cowplot)
library(magick)
library(ggh4x)
library(plotrix)
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
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE02/Primary_experiments/"
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

id_split <- function(element) {
    function(x) {
        extract <- strsplit(x, "_")[[1]][element]
        ifelse(element == 1, extract, extract <- as.numeric(extract))
        extract
    }
}

splittingSeeds <- function(value1, value2, separator) {
    value1Split <- strsplit(value1, separator)[[1]]
    value2Split <- strsplit(value2, separator)[[1]]
    #print(c(value1Split, value2Split))
    #flush.console()
    if(any(value1Split %in% value2Split)) {
        #print(value1)
        commonVal <- value1Split[value1Split %in% value2Split]
        if(length(commonVal) > 1) {return(paste(commonVal, collapse = ","))} else {
            return(commonVal)} } else {
                return(value2)}
}

marking.unresolvable.haplotype.transitions <- function(df) {
    print(strsplit(df$id[1], "_")[[1]][1])
    flush.console()
    df <- as.data.frame(df)
    byChr <- split(df, df$chr)
    corrected_chroms <- lapply(byChr, function(chrDf) {
        #print(chrDf$chr[1])
        #flush.console()
        counts <- rle(chrDf$mfh)
        if(length(counts$values) > 1){
            splitSeed <- unlist(lapply(1:length(counts$values), function(x) {
                value1Split <- strsplit(counts$values[x], ",")[[1]]
                value2Split <- strsplit(counts$values[x+1], ",")[[1]]
                if(any(value1Split %in% value2Split)) {
                    return(1)} else{ 
                        return(0)}
            } ) )
            correctSeedIdx <- unlist(lapply(1:(length(splitSeed)-1), function(x) {
                if(splitSeed[x] == 1 && splitSeed[x+1] == 0) {
                    return(x+1)} 
            } ) )
            if(length(correctSeedIdx) > 0) {splitSeed[correctSeedIdx] <- 1} else{
                splitSeed <- splitSeed}
            binSeed <- rle(splitSeed)
            summing_ties <- cumsum(binSeed$lengths)
            summing_ties2 <- c(1, cumsum(binSeed$lengths) + 1)
            first_flank <- min(which(binSeed$values == 1))
            if(first_flank == Inf) {
                #print(c(chrDf$chr, chrDf$mfh))
                #flush.console()
            }
            if(first_flank != Inf) { 
                if(first_flank == 2) {flank1 <- binSeed$lengths[1] + 1}
                if(first_flank == 1) {flank1 <- 1}
                if(flank1 != 1) {
                    flanking1 <- summing_ties2[seq(2, length(summing_ties2)-1, 2)]
                    flanking2 <- summing_ties[seq(2, length(summing_ties), 2)]  
                }
                if(flank1 == 1) {
                    flanking1 <- summing_ties2[seq(1, length(summing_ties2)-1, 2)]
                    flanking2 <- summing_ties[seq(1, length(summing_ties), 2)]  
                }
                all_flanks <- sort(c(flanking1, flanking2))
                correcting <- lapply(seq(1,length(all_flanks), 2), function(x) {
                    left_flank <- all_flanks[x]
                    right_flank <- all_flanks[x+1]
                    hapCalls <- counts$values[left_flank:right_flank]
                    tiedIdxs <- left_flank:right_flank
                    commonHaps <- Reduce(function(val1, val2) splittingSeeds(val1, val2, ","), hapCalls, accumulate = TRUE)
                    commonHaps <- commonHaps[-1]
                    #print(commonHaps)
                    #flush.console()
                    if(length(commonHaps) > 1) {
                        findTransitions <- unlist(lapply(1:(length(commonHaps)-1), function(f) {
                            val1split <- strsplit(commonHaps[f], ",")[[1]]
                            val2split <- strsplit(commonHaps[f+1], ",")[[1]]
                            if(any(val1split %in% val2split)) {
                                return(1)} else{ 
                                    return(0)}
                        } ) )
                        if(0 %in% findTransitions) {
                            cutoff <- which(findTransitions == 0) + 1
                            transitionHaps <- commonHaps[cutoff]
                            if(length(cutoff) > 1) {
                                print(chrDf$chr[1])
                                flush.console()
                                multTransitions <- unlist(lapply(1:length(cutoff), function(p) {
                                    unknownHaps <- strsplit(transitionHaps[p], ",")[[1]]
                                    findUnknowns <- unlist(lapply(hapCalls[-1], function(g) {
                                        hapsSplit <- strsplit(g, ",")[[1]]
                                        if(any(unknownHaps %in% hapsSplit)) {
                                            return(1)} else {return(0)}
                                    } ) )
                                    if(p == 1){
                                        findUnknowns[cutoff[p]:length(findUnknowns)] <- 0} else{
                                            findUnknowns[c(1:cutoff[p-1], cutoff[p]:length(findUnknowns))] <- 0
                                        }
                                    unknownIdxs <- which(findUnknowns == 1) + 1
                                    unknownIdxs <- tiedIdxs[unknownIdxs]
                                    unknownIdxs
                                } ) )
                                return(data.frame(values_idx = multTransitions, correct_founder = rep("?", length(multTransitions)), stringsAsFactors = FALSE))
                            } else {
                                unknownHaps <- strsplit(transitionHaps, ",")[[1]]
                                findUnknowns <- unlist(lapply(hapCalls[-1], function(g) {
                                    hapsSplit <- strsplit(g, ",")[[1]]
                                    if(any(unknownHaps %in% hapsSplit)) {
                                        return(1)} else {return(0)}
                                } ) )
                                findUnknowns[cutoff:length(findUnknowns)] <- 0
                                unknownIdxs <- which(findUnknowns == 1) + 1
                                unknownIdxs <- tiedIdxs[unknownIdxs]
                                return(data.frame(values_idx = unknownIdxs, correct_founder = rep("?", length(unknownIdxs)), stringsAsFactors = FALSE)) } 
                        }
                    }
                } )
                correcting1 <- Filter(Negate(function(i) is.null(unlist(i))), correcting)
                correcting_df <- do.call(rbind, correcting1)
                counts$values[correcting_df$values_idx] <- correcting_df$correct_founder
                corrected_mfh <- unlist(lapply(1:length(counts$lengths), function(x) {
                    rep(counts$values[x], counts$lengths[x])
                } ) )
                chrDf$mfh <- corrected_mfh
            } else {
                chrDf$mfh <- chrDf$mfh } } else {
                    chrDf$mfh <- chrDf$mfh}
        chrDf
    } )
    correctedData <- do.call(rbind, corrected_chroms)
}

adjust.multiple.haplotype.calls <- function(unknownDf) {
    print(strsplit(unknownDf$id[1], "_")[[1]][1])
    flush.console()
    byChr <- split(unknownDf, unknownDf$chr)
    corrected_chroms <- lapply(byChr, function(chrDf) {
        print(chrDf$chr[1])
        flush.console()
        counts <- rle(chrDf$mfh)
        if(length(counts$values) > 1){
            splitSeed <- unlist(lapply(1:length(counts$values), function(x) {
                value1Split <- strsplit(counts$values[x], ",")[[1]]
                value2Split <- strsplit(counts$values[x+1], ",")[[1]]
                if(any(value1Split %in% value2Split)) {
                    return(1)} else{ 
                        return(0)}
            } ) )
            correctSeedIdx <- unlist(lapply(1:(length(splitSeed)-1), function(x) {
                if(splitSeed[x] == 1 && splitSeed[x+1] == 0) {
                    return(x+1)} 
            } ) )
            if(length(correctSeedIdx) > 0) {
                fixOne <- unlist(lapply(1:length(correctSeedIdx), function(x) {
                    if(correctSeedIdx[x] < length(splitSeed)) {
                        if(splitSeed[correctSeedIdx[x] + 1] == 1) { 
                            return(correctSeedIdx[x] + 1)
                        }
                    }
                } ) )
                splitSeed[correctSeedIdx] <- 1
                splitSeed[fixOne] <- 2
            } else{
                splitSeed <- splitSeed}
            splitSeedDT <- data.table(seed = splitSeed)[, row := .I][, leaveOut := seed != 2][, seedGrp := rleid(seed)*leaveOut]
            findGrp <- splitSeedDT[seedGrp == 0, row] + 1
            findGrp2 <- findGrp[findGrp <= splitSeedDT[, .N]]
            splitSeedDT2 <- splitSeedDT[, seedMult := ifelse(seed == 0, 0, 1)]
            splitSeedDT2 <- splitSeedDT[seedGrp == 0, seedGrp := splitSeedDT[findGrp2, seedGrp]][, seedGrp := seedGrp*seedMult]
            splitSeedDT2[, haps := counts$values][, hapLengths := counts$lengths][, numCommas := lengths(regmatches(haps, gregexpr(",", haps)))]
            grpCollapseDT <- splitSeedDT2[seedGrp != 0, .(collapsedLengths = sum(hapLengths)), by = seedGrp]
            minHaps <- splitSeedDT2[seedGrp != 0, .(minHaps = haps[numCommas == min(numCommas)]), by = seedGrp]
            commonHaps <- minHaps[, Reduce(function(val1, val2) splittingSeeds(val1, val2, ","), minHaps), by = seedGrp]
            hapsResolved <- merge(grpCollapseDT, commonHaps, by = "seedGrp")
            if(hapsResolved[, .N] > 0) {
                setnames(hapsResolved, old = c("collapsedLengths", "V1"), new = c("hapLengths", "haps"))
            }
            findTiedRows <- splitSeedDT2[seedGrp != 0,.(row = row[1]), by = seedGrp]
            hapsResolved[, row := findTiedRows$row]
            newSeedDT <- splitSeedDT2[seedGrp == 0][, c("row", "seedGrp", "haps", "hapLengths")]
            if(hapsResolved[, .N] > 0) {
                correctedSeedDT <- rbind(newSeedDT, hapsResolved)[order(row)]
            } else {
                correctedSeedDT <- newSeedDT
            }
            correctedMFH <- rep(correctedSeedDT$haps, correctedSeedDT$hapLengths)
            chrDf$mfh <- correctedMFH
        } else {
            chrDf$mfh <- chrDf$mfh}
        chrDf
    } )
    correctedData <- do.call(rbind, corrected_chroms)
    correctedData
}

infer.haplotype.calls <- function(adjustedDf) {
    print(strsplit(adjustedDf$id[1], "_")[[1]][1])
    flush.console()
    byChr <- split(adjustedDf, adjustedDf$chr)
    corrected_chroms <- lapply(byChr, function(chrDf) {
        print(chrDf$chr[1])
        flush.console()
        counts <- rle(chrDf$mfh)
        orig_counts <- counts
        unknowns_in_counts <- grep('\\?', counts$values)
        if(length(counts$values) > 2){
            ties_df <- data.frame(ties = rep(0, length(counts$values)), stringsAsFactors = FALSE)
            ties_df$ties[unknowns_in_counts] <- 1
            bin_ties <- rle(ties_df$ties)
            summing_ties <- cumsum(bin_ties$lengths)
            summing_ties2 <- cumsum(bin_ties$lengths) + 1
            first_flank <- min(which(bin_ties$values == 0))
            last_flank <- max(which(bin_ties$values == 1))
            if(first_flank != Inf & last_flank > first_flank) { 
                ##########  Error here fix #################
                flanking1 <- summing_ties[seq(first_flank, last_flank, 2)]
                flanking2 <- summing_ties2[seq(first_flank + 1, last_flank + 1, 2)]
                all_flanks <- sort(c(flanking1, flanking2))
                counts$lengths <- cumsum(counts$lengths)
                correcting <- lapply(seq(1,length(all_flanks), 2), function(x) {
                    left_flank <- all_flanks[x]
                    right_flank <- all_flanks[x+1]
                    left_flank_id <- counts$values[left_flank]
                    right_flank_id <- counts$values[right_flank]
                    split_left <- strsplit(left_flank_id, ",")[[1]]
                    split_right <- strsplit(right_flank_id, ",")[[1]]
                    unknown <- (left_flank + 1):(right_flank - 1)
                    if(any(split_left %in% split_right)){
                        inferredCall <- split_left[split_left %in% split_right]
                        if(length(inferredCall > 1)) {
                            inferredCall <- paste(inferredCall, collapse = ",")
                        } else {
                            inferredCall <- inferredCall
                        }
                        return(data.frame(values_idx = unknown, correct_founder = rep(inferredCall, length(unknown)), stringsAsFactors = FALSE))}
                } )
                correcting1 <- Filter(Negate(function(i) is.null(unlist(i))), correcting)
                correcting_df <- do.call(rbind, correcting1)
                orig_counts$values[correcting_df$values_idx] <- correcting_df$correct_founder
                corrected_mfh <- unlist(lapply(1:length(orig_counts$lengths), function(x) {
                    rep(orig_counts$values[x], orig_counts$lengths[x])
                } ) )
                chrDf$mfh <- corrected_mfh
            } else {
                chrDf$mfh <- chrDf$mfh} } else {
                    chrDf$mfh <- chrDf$mfh}
        chrDf
    } )
    correctedData <- do.call(rbind, corrected_chroms)
}


defining.haplotype.blocks <- function(inferred_hap_df) {
    print(strsplit(inferred_hap_df$id[1], "_")[[1]][1])
    flush.console()
    byChr <- split(inferred_hap_df, inferred_hap_df$chr)
    corrected_chroms <- lapply(byChr, function(chrDf) {
        print(chrDf$chr[1])
        flush.console()
        calls <- rle(chrDf$mfh)
        startInterval <- c(1, cumsum(calls$lengths[-length(calls$lengths)]) + 1)
        endInterval <- cumsum(calls$lengths)
        call_sums <- sort(c(startInterval, endInterval))
        intervals <- chrDf[call_sums,]
        intervals
    } )
    correctedData <- do.call(rbind, corrected_chroms)
}

sizing.haplotype.blocks <- function(hapIntervalDf) {
    print(strsplit(hapIntervalDf$id[1], "_")[[1]][1])
    flush.console()
    hap_blck_lengths <- unlist(lapply(seq(1,nrow(hapIntervalDf),2), function(x) {
        #print(x)
        if(!grepl('\\?', hapIntervalDf$mfh[x])) {
            if(!is.na(hapIntervalDf$chr[x+2])) {
                if(hapIntervalDf$chr[x+2] == hapIntervalDf$chr[x]) {
                    distance <- hapIntervalDf$pos[x+2] - hapIntervalDf$pos[x]
                    if(distance < 0) {print(x)}
                    return(distance)} else {
                        distance <- hapIntervalDf$pos[x+1] + 1000 - hapIntervalDf$pos[x]
                        if(distance < 0) {print(x)}
                        return(distance)}
            } else if(is.na(hapIntervalDf$chr[x+2])) {
                distance <- hapIntervalDf$pos[x+1] + 1000 - hapIntervalDf$pos[x]
                return(distance)} 
        }
    } ) )
    hap_blck_lengths
}

sizing.unknown.blocks <- function(hapIntervalDf) {
    print(strsplit(hapIntervalDf$id[1], "_")[[1]][1])
    flush.console()
    hap_blck_lengths <- unlist(lapply(seq(1,nrow(hapIntervalDf),2), function(x) {
        #print(x)
        if(grepl('\\?', hapIntervalDf$mfh[x])) {
            if(!is.na(hapIntervalDf$chr[x+2])) {
                if(hapIntervalDf$chr[x+2] == hapIntervalDf$chr[x]) {
                    distance <- hapIntervalDf$pos[x+2] - hapIntervalDf$pos[x]
                    if(distance < 0) {print(x)}
                    return(distance)} else {
                        distance <- hapIntervalDf$pos[x+1] + 1000 - hapIntervalDf$pos[x]
                        if(distance < 0) {print(x)}
                        return(distance)}
            } else if(is.na(hapIntervalDf$chr[x+2])) {
                distance <- hapIntervalDf$pos[x+1] + 1000 - hapIntervalDf$pos[x]
                return(distance)} 
        }
    } ) )
    sum(hap_blck_lengths)
}

sizing.mult.blocks <- function(hapIntervalDf) {
    print(strsplit(hapIntervalDf$id[1], "_")[[1]][1])
    flush.console()
    hap_blck_lengths <- unlist(lapply(seq(1,nrow(hapIntervalDf),2), function(x) {
        #print(x)
        if(grepl(',', hapIntervalDf$mfh[x])) {
            if(!is.na(hapIntervalDf$chr[x+2])) {
                if(hapIntervalDf$chr[x+2] == hapIntervalDf$chr[x]) {
                    distance <- hapIntervalDf$pos[x+2] - hapIntervalDf$pos[x]
                    if(distance < 0) {print(x)}
                    return(distance)} else {
                        distance <- hapIntervalDf$pos[x+1] + 1000 - hapIntervalDf$pos[x]
                        if(distance < 0) {print(x)}
                        return(distance)}
            } else if(is.na(hapIntervalDf$chr[x+2])) {
                distance <- hapIntervalDf$pos[x+1] + 1000 - hapIntervalDf$pos[x]
                return(distance)} 
        }
    } ) )
    hap_blck_lengths
}

number.haplotype.blocks <- function(hapIntervalDf) {
    print(strsplit(hapIntervalDf$id[1], "_")[[1]][1])
    flush.console()
    hapBlcks <- hapIntervalDf[-(grep('\\?', hapIntervalDf$mfh)),]
    nrow(hapBlcks)/2
}

detecting.recombination.events <- function(hapIntervalDf) {
    print(strsplit(hapIntervalDf$id[1], "_")[[1]][1])
    flush.console()
    hapBlcks <- hapIntervalDf[-(grep('\\?', hapIntervalDf$mfh)),]
    byChr <- split(hapBlcks, hapBlcks$chr)
    detectingRecEvents <- unlist(lapply(byChr, function(chrDf) {
        print(chrDf$chr[1])
        flush.console()
        if(length(unique(chrDf$mfh)) > 1){
            events <- unlist(lapply(seq(1, nrow(chrDf)-2, 2), function(x) {
                hap1 <- strsplit(chrDf$mfh[x], ",")[[1]]
                hap2 <- strsplit(chrDf$mfh[x+2], ",")[[1]]
                if(any(hap1 %in% hap2)) {return(0)} else {
                    return(1)}
            } ) )
        }
    } ) )
    length(detectingRecEvents[detectingRecEvents == 1])
}

finding.haplotypes.present <- function(hapIntervalDf) {
    hapsPresent <- unique(hapIntervalDf$mfh)
    unknowns <- grep('\\?', hapsPresent)
    multipleFounders <- grep(',', hapsPresent)	
    hapsPresent <- hapsPresent[-c(unknowns, multipleFounders)]
    length(hapsPresent)
}

finding.haplotypes.present.ind.chr <- function(hapIntervalDf, chrom) {
    hapIntervalDf <- hapIntervalDf[hapIntervalDf$chr == chrom,]
    hapsPresent <- unique(hapIntervalDf$mfh)
    unknowns <- grep('\\?', hapsPresent)
    multipleFounders <- grep(',', hapsPresent)	
    hapsPresent <- hapsPresent[-c(unknowns, multipleFounders)]
    #if("AB1" %in% hapsPresent) {print(hapIntervalDf$pos[hapIntervalDf$mfh == "AB1"])}
    #if("A8" %in% hapsPresent) {print(c("A8",hapIntervalDf$id[1]))}
    hapsPresent
}


plotting.haploid.haplotypes <- function(rhc_interval_dfs, fname) {
    offsets <- read.table("newoffsets.txt", header = T)
    offsets$lines <- ceiling(offsets[,2]/50)
    offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
    offsets$chr <- 1:17
    offsets <- offsets[-nrow(offsets),] # no mitochondira
    g_l <- c(0, cumsum(offsets$len))
    ch.bounds <- c(0, g_l[1:16] + offsets[,2])
    #correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
    max.pos <- g_l[length(g_l)]
    mid.ch <- diff(ch.bounds)/2
    midpt.ch <- ch.bounds[2:17] - mid.ch
    plotting_factors <- c(topHaps, 'unknown', 'unresolvable')
    plotting_colors <- c(mycols, 'darkgray')
    color_df <- data.frame(factors = plotting_factors, colors = plotting_colors, stringsAsFactors = FALSE)
    title <- strsplit(fname, "/")[[1]]
    pTitle <- title[[length(title)]]
    week1 <- rhc_interval_dfs[[1]]$week[1]
    week2 <- rhc_interval_dfs[[length(rhc_interval_dfs)]]$week[1]
    weekChange <- unlist(lapply(rhc_interval_dfs, function(x) {
        if(x$week[1] == week2 & week1 != week2) {return(x$id[1])}
    }))
    if(length(weekChange > 0)) {
        changeIdx <- weekChange[1]} else {
            changeIdx <- "nothing"
        }
    pdf(paste0(fname, '.pdf'), height = 4, width = 6)
    par(mar = c(0,0,0,6), oma = c(3,2,3,2), xpd = TRUE)
    par(mgp = c(1.7,0.75,0))
    if(week1 == week2) {
        layout(matrix(1:length(rhc_interval_dfs), length(rhc_interval_dfs), byrow = TRUE)) } else {
            layout(matrix(1:(length(rhc_interval_dfs)+1), length(rhc_interval_dfs)+1, byrow = TRUE)) }
    counts <- 0
    plot_all <- lapply(rhc_interval_dfs, function(df) {
        counts <<- counts + 1
        if(df$Chemical[1] == "base") {
            id <- substr(df$id[1], nchar(df$id[1])-2, nchar(df$id[1]))} else {
                id <- paste0("wk", substr(df$id[1], 4, 5), substr(df$id[1], nchar(df$id[1])-2, nchar(df$id[1])))
            }
        print(id)
        flush.console()
        df$mfh[df$mfh == '?'] <- 'unresolvable' 
        #df$mfhp <- df$mfh
        df$mfh <- gsub(",", "", df$mfh)
        df$mfh[df$mfh != 'unresolvable' & !df$mfh %in% color_df$factors] <- 'unknown'
        df$cols <- unlist(lapply(df$mfh, function(x) {
            color_df$colors[color_df$factors == x]
        } ) )
        #par(mar = c(1.25,0,0,0), mgp = c(1.2,1,0), oma = c(2,4,3,1))
        #par(mar = c(0,0,0,0), oma = c(1,0.25,0,0.5), mgp = c(0.75,0.5,0))
        xmin <- 0
        xmax <- max.pos
        plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(xmin, xmax), ylim= c(0, 1), axes="FALSE", yaxs="i", xaxs="i")
        if(df$id[1] != changeIdx) {
            plotting <- lapply(seq(1,nrow(df), 2), function(x){
                rect(df$gp[x], 0.1, df$gp[x+1], 0.9, density = NA, col = df$cols[x], border = NA)
            } )
            axis(2, at = 0.5, labels = id, las = 1, cex.axis = 0.5, lwd = 0, line = -0.5)
            clip(0,xmax,0,1)
            #abline(h= 0, col = 'black', lwd = 1)
            abline(v= ch.bounds, col = 'black', lwd = 0.5, lty = 1)
        }
        if(counts == 1) {
            abline(h = 1, col = 'black', lwd = 1)
            #title(main = pTitle, cex.main = 1, line = 0.5, xpd = NA)
            legend(grconvertX(1.005, "npc", "user"),grconvertY(1.25, "npc", "user"), legend = unique(plotting_factors), col = plotting_colors, lwd = 2, cex = 0.75, xpd = NA, bty = "n") }
        if(counts == length(rhc_interval_dfs)) {
            axis(1, at = midpt.ch, labels = as.roman(1:16), las = 1, cex.axis = 0.75, lwd = 0, line = -0.5)
            clip(0,xmax,0,1)
            abline(h= 0, col = 'black', lwd = 1)}
    } )		
    dev.off()
}

plotting.zoomed.in.haploid.haplotypes <- function(rhc_interval_dfs, fname, chroms) {
    offsets <- read.table("newoffsets.txt", header = T)
    offsets$lines <- ceiling(offsets[,2]/50)
    offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
    offsets$chr <- 1:17
    g_l <- c(0, cumsum(offsets$len[chroms]))
    ch.bounds <- c(0, g_l[1:length(chroms)] + offsets[chroms,2])
    #correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
    max.pos <- g_l[length(g_l)]
    mid.ch <- diff(ch.bounds)/2
    midpt.ch <- ch.bounds[2:length(ch.bounds)] - mid.ch
    colPal <- c("dodgerblue2", "powderblue", # red
                "lemonchiffon2",
                "#6A3D9A", # purple
                "#FF7F00", # orange
                "gold1",
                "skyblue2", "#FB9A99", # lt pink
                "palegreen2",
                "#CAB2D6", # lt purple
                "tan1",
                "coral2",
                "tomato", # lt orange
                "maroon", "blue1", "steelblue4", "darkturquoise")
    plotting_factors <- c('unknown', 'multiple_founders', newFounderNames)
    plotting_colors <- c('darkgray', 'black', colPal)
    color_df <- data.frame(factors = plotting_factors, colors = plotting_colors, stringsAsFactors = FALSE)
    #pdf(paste0(fname, '.pdf'), height = 4, width = 6)
    #dev.new(height = 4, width = 6)
    pdf(paste0(fname, '.pdf'), height = 4, width = 6)
    par(mar = c(0,0,0,6), oma = c(5,2,2,2), xpd = TRUE)
    par(mgp = c(1.7,0.75,0))
    layout(matrix(1:length(rhc_interval_dfs), length(rhc_interval_dfs), byrow = TRUE))
    counts <- 0
    plot_all <- lapply(rhc_interval_dfs, function(df) {
        counts <<- counts + 1
        id <- paste0('rhc ', counts)
        df <- df[df$chr == chroms,]
        df$gp <- df$pos
        df$mfh[df$mfh == '?'] <- 'unknown' 
        df$mfhp <- df$mfh
        df$mfhp[grepl(',', df$mfhp)] <- 'multiple_founders'
        df$cols <- unlist(lapply(df$mfhp, function(x) {
            color_df$colors[color_df$factors == x]
        } ) )
        #par(mar = c(1.25,0,0,0), mgp = c(1.2,1,0), oma = c(2,4,3,1))
        #par(mar = c(0,0,0,0), oma = c(1,0.25,0,0.5), mgp = c(0.75,0.5,0))
        xmin <- 0
        xmax <- max.pos
        plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(xmin, xmax), ylim= c(0, 1), axes="FALSE", yaxs="i", xaxs="i", bty = 'n')
        plotting <- lapply(seq(1,nrow(df), 2), function(x){
            rect(df$gp[x], 0.25, df$gp[x+1], 0.75, density = NA, col = df$cols[x], border = NA)
        } )
        axis(2, at = 0.5, labels = id, las = 1, cex.axis = 0.75, lwd = 0, line = -0.5)
        clip(0,xmax,0,1)
        abline(h= 0, col = 'black', lwd = 1)
        abline(v= ch.bounds, col = 'black', lwd = 1, lty = 1)
        if(counts == 1) {
            abline(h= 1, col = 'black', lwd = 1)
            title(main = fname, cex.main = 1, line = 0.5, xpd = NA)
            legend(grconvertX(1.005, "npc", "user"), grconvertY(1.25, "npc", "user"), legend = unique(plotting_factors), col = plotting_colors, lwd = 2, cex = 0.75, xpd = NA, bty = "n") }
        if(counts == length(rhc_interval_dfs)) { 
            axis(1, at = seq(0, df$gp[nrow(df)], 100000), labels = seq(0, df$gp[nrow(df)]/1000, 100), las = 1, cex.axis = 0.75, line = 1, lwd = 0.75)
            #par(mar = c(3,0,0,6))
            mtext("Position (kb)", side = 1, line = 3, cex = 0.5, xpd = NA)
        }
    } )		
    dev.off()
}

plot.haplotype.sizes <- function(list_of_hap_sizes, sample_names) {
    png("Haplotype_sizes_histograms.png", height = 6, width = 6, units = "in", res = 150)
    #dev.new(height = 8, width = 8)
    par(mar = c(3,4.5,2,0.5), mgp = c(1.2,1,0), oma = c(0,0,0,0))
    par(mfrow=c(7,1))
    counter <- 0
    plotting <- lapply(list_of_hap_sizes, function(x) {
        counter <<- counter + 1
        sampleHist <- hist(x, breaks = 25, plot = FALSE)
        sampleHist$density <- sampleHist$counts/sum(sampleHist$counts)*100
        breakDivide <- 1000
        sampleHist$breaks <- sampleHist$breaks/breakDivide
        print(sampleHist$breaks)
        flush.console()
        plot(sampleHist, freq = FALSE, main = "", xlim = c(0, 1500), ylim = c(0,50), xlab = "", ylab = "", las = 1)
        mtext("Haplotype block length (kb)", side = 1, line = 2, cex = 0.75)
        mtext("Percent of haplotype blocks in bin", side = 2, line = 2.5, cex = 0.75)
        text(grconvertX(-0.16, "npc", "user"),grconvertY(1.07, "npc", "user"), bquote(bold(.(LETTERS[counter]))), family = "sans", cex = 1.5, xpd = NA)
    } )
    dev.off()
}

plot.unknown.sizes <- function(list_of_unknown_sizes, sample_names) {
    pdf("Unknown_haplotype_sizes_histograms.pdf", height = 8, width = 8)
    #dev.new(height = 8, width = 8)
    par(mar = c(4,5,2,2), mgp = c(1.2,1,0), oma = c(0,0,0,0))
    par(mfrow=c(2,1))
    counter <- 0
    plotting <- lapply(list_of_unknown_sizes, function(x) {
        counter <<- counter + 1
        sampleHist <- hist(x, breaks = seq(0,200000,5000), plot = FALSE)
        sampleHist$density <- sampleHist$counts/sum(sampleHist$counts)*100
        breakDivide <- 1000
        sampleHist$breaks <- sampleHist$breaks/breakDivide
        print(sampleHist$breaks)
        flush.console()
        plot(sampleHist, freq = FALSE, main = paste0(sample_names[counter], 'Unknown haplotype block size distribution'), xlim = c(0, 200), ylim = c(0,60), xlab = "", ylab = "")
        mtext("Unknown haplotype block length (kb)", side = 1, line = 2, cex = 1)
        mtext("Percent of haplotype blocks in bin", side = 2, line = 2.5, cex = 1)
        text(grconvertX(-0.085, "npc", "user"),grconvertY(1.07, "npc", "user"), bquote(bold(.(LETTERS[counter]))), family = "sans", cex = 1.5, xpd = NA)
    } )
    dev.off()
}

# ============================================================================
# Load data needed for downstream analyses

founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)

#indHapDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Individual_tx_tables/", analysisType = "hap_freqs_250snps", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]|^BAS02")
indHapDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Individual_tx_tables/", analysisType = "hap_freqs", samplePattern = "^SEE[0-9][0-9]B02[A-Z][A-Z][0-9][0-9][0-9]R[0-9][0-9]|^BAS02")

allBind <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_All_Reps_Hets.txt"))
topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))
#hetDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"), header = T)
treatKeyDT <- fread("treatment_key.txt", header = T)
offsets <- read.table("newoffsets.txt", header = T)

offsets$lines <- ceiling(offsets[,2]/50)
offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
offsets$chr <- 1:16
g_l <- c(0, cumsum(offsets$len))
ch.bounds <- c(0, g_l[1:16] + offsets[,2])
#correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
max.pos <- g_l[length(g_l)]
mid.ch <- diff(ch.bounds)/2
midpt.ch <- ch.bounds[2:17] - mid.ch
gray <- col2rgb('grey50')
chromBarCol <- rgb(gray[1],gray[2],gray[3], maxColorValue = 255, alpha = 100)
colorCodes = c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
hx = sapply(strsplit(colorCodes, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
hx2 = hx[-c(5,21,24)]   
mycols = c(hx2,"#000000")
mycolsNames <- unlist(lapply(mycols,  function(x) color.id(x)[1]))
mycols2 <- c("blue", "red")
myColPal <- add.alpha(mycols2, 0.5)

# ============================================================================
# Data wrangle frequency change of all haplotypes (for each replicate), with each replicate a different shape and the haplotypes colored consistently with previous plots. 

### Looks like some 'haploids' are diploids- check raw data!

popDT <- do.call(rbind, c(indHapDTs, list(fill = TRUE)))
popDT <- popDT[grepl("H[0-9][0-9]$", id)]
chemSplit <- split(popDT, popDT$id)

chemLooper <- lapply(chemSplit, function(chem) {
    print(chem$id[1])
    flush.console()
    freqs <- chem[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
    collHaps <- chem[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
    allHaps <- cbind(chem$chr, chem$pos, chem$gp, chem$Chemical, chem$Replicate, chem$id, collHaps, freqs)
    setnames(allHaps, c("chr", "pos", "gp", "Chemical", "Replicate", "id", paste0("V", 1:34)))
    freqsDF <- as.data.frame(freqs)
    newFreqs <- lapply(freqsDF, function(y) gsub("^-\\d.*|^\\d.*", 1, y))
    newFreqsDF <- as.data.frame(do.call(cbind, newFreqs))
    newFreqsDF2 <- sapply(newFreqsDF, as.numeric)
    allHaps$reps <- rowSums(newFreqsDF2, na.rm = T)
    chroms <- rep(unlist(allHaps$chr), allHaps$reps)
    positions <- rep(unlist(allHaps$pos), allHaps$reps) 
    gpositions <- rep(unlist(allHaps$gp), allHaps$reps)
    replicates <- sort(rep(unlist(allHaps$Replicate), allHaps$reps))
    allfreqs <- as.numeric(unlist(t(freqs)))
    allfreqs2 <- na.omit(allfreqs)
    hapVec <- as.character(unlist(t(collHaps)))
    hapVec2 <- na.omit(hapVec)
    idxDT <- data.table(chr = chroms, pos = positions, gp = gpositions, Chemical = allHaps$Chemical[1], Replicate = replicates, id = allHaps$id[1], haps = hapVec2, freqs = allfreqs2)
    idxDT[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
    #hapColors <- match(idxDT$haps, topHaps)
    #hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
    #Colors <- mycols[hapColors]
    #Colors2 <- mycolsNames[hapColors]
    #idxDT[, "color.codes" := .(Colors)]
    #idxDT[, "color.names" := .(Colors2)]
    mfhDT <- idxDT[idxDT[, .I[freqs == max(freqs)], by=gp]$V1]
    mfhDT$mfh <- as.character(mfhDT$haps)
    mfhDT$mfh <- gsub('([0-9])([A-B])', '\\1,\\2', mfhDT$haps)
    #hetDT <- mfhDT[freqs >= 0.5 & freqs <= 0.6] # Keeping potentially aneuploid regions
    mfhDT <- mfhDT[freqs < 0.9, mfh := '?'] # keeping all regions greater than 0.9 (should be haploid)
    #mfhDT[gp %in% hetDT$gp, mfh := hetDT$mfh]
    #chem3 <- chem2[freqs >= 0.95]
    #colorMappings <- chem3[, unique(haps)]
    mfhDT
} )



##############  Inferring tied haplotype calls  ####################

markUnknownTransitions <- Map(marking.unresolvable.haplotype.transitions, chemLooper)
adjustedHapCalls <- Map(adjust.multiple.haplotype.calls, markUnknownTransitions)
inferredHapCalls <- Map(infer.haplotype.calls, adjustedHapCalls)

####################################  Delineating haplotype blocks in RHCs  ###########################################

definedHapBlocks <- Map(defining.haplotype.blocks, inferredHapCalls)

####################################  Number of haplotype blocks in RHCs  ###########################################

allHapBlocks <- do.call(rbind, definedHapBlocks)
allHapBlocks$week <- substr(allHapBlocks$id, 4,5)
allHapBlocks$week[allHapBlocks$Chemical == "base"] <- 0
allHapBlocks$chemWeek <- paste0(allHapBlocks$Chemical, "_wk", allHapBlocks$week)
allHapBlocks$chemWeek <- gsub(" ", "_", allHapBlocks$chemWeek)
allHapsSplit <- split(allHapBlocks, allHapBlocks$id)

hapLengths <- lapply(allHapsSplit, function(chem) {
    #fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/", chem$Chemical[1], "_50snp_window_filtered")
    unresolvLength <- sizing.unknown.blocks(chem)
    fracUnresolv <- unresolvLength/max.pos
    if(fracUnresolv < 0.25) {
        numHapBlocks <- number.haplotype.blocks(chem)
        numRecEvents <- detecting.recombination.events(chem)
        return(data.table(id = chem$id[1], hapBlocks = numHapBlocks, recEvents = numRecEvents, chemical = chem$Chemical[1], week = chem$week[1], chemWeek = chem$chemWeek[1]))
    }
} )

hapBlocks <- do.call(rbind, hapLengths)
avgMedBlocks <- hapBlocks[, .(avgBlocks = mean(hapBlocks), medBlocks = median(hapBlocks)), by = chemWeek]

####################################  Size and number  of haplotype blocks in RHCs  ###########################################

allHapBlocks <- do.call(rbind, definedHapBlocks)
allHapBlocks$week <- substr(allHapBlocks$id, 4,5)
allHapBlocks$week[allHapBlocks$Chemical == "base"] <- 0
allHapBlocks$chemWeek <- paste0(allHapBlocks$Chemical, "_wk", allHapBlocks$week)
allHapBlocks$chemWeek <- gsub(" ", "_", allHapBlocks$chemWeek)
allHapsSplit <- split(allHapBlocks, allHapBlocks$id)

hapLengths <- lapply(allHapsSplit, function(chem) {
    #fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/", chem$Chemical[1], "_50snp_window_filtered")
    unresolvLength <- sizing.unknown.blocks(chem)
    fracUnresolv <- unresolvLength/max.pos
    if(fracUnresolv < 0.25) {
        numHapBlocks <- number.haplotype.blocks(chem)
        numRecEvents <- detecting.recombination.events(chem)
        hapBlockLengths <- sizing.haplotype.blocks(chem)
        return(data.table(id = chem$id[1], hapBlocks = numHapBlocks, recEvents = numRecEvents, hapLengths = hapBlockLengths, chemical = chem$Chemical[1], week = chem$week[1], chemWeek = chem$chemWeek[1]))
    }
} )

hapBlocks <- do.call(rbind, hapLengths)
wk11medBlcks <- hapBlocks[grepl("wk11", chemWeek), median(hapLengths)]
avgMedBlocks <- hapBlocks[, .(avgBlocks = mean(hapBlocks), medBlocks = median(hapBlocks), avgBlcLens = mean(hapLengths), medBlcLens = median(hapLengths), avgRecEvents = mean(recEvents), medRecEvents = median(recEvents)), by = chemWeek]

allIDs <- hapBlocks[, .(clones = uniqueN(id)), by = chemWeek]
allDT <- avgMedBlocks[allIDs, on = "chemWeek"]

###########################  Plotting haplotype block size distributions  ############

plot.haplotype.sizes <- function(list_of_hap_sizes, sample_names) {
    fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/Haplotype_sizes_histograms.png")
    png(fileName, height = 6, width = 6, units = "in", res = 150)
    #dev.new(height = 8, width = 8)
    par(mar = c(3,4.5,2,0.5), mgp = c(1.2,1,0), oma = c(0,0,0,0))
    par(mfrow=c(7,1))
    counter <- 0
    plotting <- lapply(list_of_hap_sizes, function(dt) {
        counter <<- counter + 1
        x <- dt$hapLengths
        sampleHist <- hist(x, breaks = 25, plot = FALSE)
        sampleHist$density <- sampleHist$counts/sum(sampleHist$counts)*100
        breakDivide <- 1000
        sampleHist$breaks <- sampleHist$breaks/breakDivide
        print(sampleHist$breaks)
        flush.console()
        plot(sampleHist, freq = FALSE, main = "", xlim = c(0, 1500), ylim = c(0,50), xlab = "", ylab = "", las = 1)
        mtext("Haplotype block length (kb)", side = 1, line = 2, cex = 0.75)
        mtext("Percent of haplotype blocks in bin", side = 2, line = 2.5, cex = 0.75)
        text(grconvertX(-0.16, "npc", "user"),grconvertY(1.07, "npc", "user"), bquote(bold(.(LETTERS[counter]))), family = "sans", cex = 1.5, xpd = NA)
    } )
    dev.off()
}

list_of_hap_sizes <- split(hapBlocks, hapBlocks$chemWeek)
sample_names <- unique(hapBlocks$chemWeek)
plottingHapSizes <- plot.haplotype.sizes(list_of_hap_sizes, sample_names)

fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/Haploid_haplotype_sizes_histograms.pdf")
panelPlot <- ggplot() + geom_histogram(data=hapBlocks, aes(x = hapLengths), binwidth = 2000) + facet_wrap(~chemWeek, ncol = 1) + theme_bw(base_size = 12) + theme(strip.text.x = element_text(color = "black", face = "plain", size = 10, margin = margin())) + theme(strip.background = element_rect(color="black", fill="gray", size=0.5, linetype="solid")) + theme(axis.text.x = element_text(angle = -45)) + coord_cartesian(xlim = c(0, 250000)) 
ggsave(file = fileName, panelPlot, width = 8.5, height = 9, units = "in")


fileName2 <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/Haploid_haplotype_sizes_boxplot_no_base.pdf")
hapPlot <- ggplot(hapBlocks[chemical != "base"], aes(x=chemWeek, y=hapLengths/1000)) + geom_boxplot(notch=FALSE) + theme_bw(base_size = 10)  + xlab("") + ylab("unrecombined block size (kb)") + coord_cartesian(ylim = c(0, 100)) + theme(axis.text.x = element_text(angle = -45))
ggsave(file = fileName2, hapPlot, width = 8.5, height = 9, units = "in")

#################################  Plotting RHC haplotypes each chemical and week separately ######################################

allHapBlocks <- do.call(rbind, definedHapBlocks)
allHapBlocks$week <- substr(allHapBlocks$id, 4,5)
allHapBlocks$week[allHapBlocks$Chemical == "base"] <- 0
allHapBlocks$chemWeek <- paste0(allHapBlocks$Chemical, "_wk", allHapBlocks$week)
allHapBlocks$chemWeek <- gsub(" ", "_", allHapBlocks$chemWeek)
allHapsSplit <- split(allHapBlocks, allHapBlocks$chemWeek)

haploidPlotter <- lapply(allHapsSplit, function(chem) {
    fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/", chem$chemWeek[1], "_250snp_window")
    chemDFs <- split(chem, chem$id)
    hapPlot <- plotting.haploid.haplotypes(chemDFs, fileName)
} )

#################################  Plotting RHC haplotypes each chemical separately, weeks together ######################################

allHapBlocks <- do.call(rbind, definedHapBlocks)
allHapBlocks$week <- substr(allHapBlocks$id, 4,5)
allHapBlocks$week[allHapBlocks$Chemical == "base"] <- 0
allHapBlocks$chemWeek <- paste0(allHapBlocks$Chemical, "_wk", allHapBlocks$week)
allHapBlocks$chemWeek <- gsub(" ", "_", allHapBlocks$chemWeek)
allHapsSplit <- split(allHapBlocks, allHapBlocks$Chemical)

haploidPlotter <- lapply(allHapsSplit, function(chem) {
    fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/", chem$Chemical[1], "_50snp_window_filtered_notitle")
    chemDFs <- split(chem, chem$id)
    filteredDFs <- lapply(chemDFs, function(x) {
        unresolvLength <- sizing.unknown.blocks(x)
        fracUnresolv <- unresolvLength/max.pos
        if(fracUnresolv < 0.25) {return(x)}
    })
    filteredDF <- do.call(rbind, filteredDFs)
    splitDFs <- split(filteredDF, filteredDF$id)
    hapPlot <- plotting.haploid.haplotypes(splitDFs, fileName)
} )

#################################  Plotting RHC haplotypes each chemical separately, weeks together, unfiltered ######################################

allHapBlocks <- do.call(rbind, definedHapBlocks)
allHapBlocks$week <- substr(allHapBlocks$id, 4,5)
allHapBlocks$week[allHapBlocks$Chemical == "base"] <- 0
allHapBlocks$chemWeek <- paste0(allHapBlocks$Chemical, "_wk", allHapBlocks$week)
allHapBlocks$chemWeek <- gsub(" ", "_", allHapBlocks$chemWeek)
allHapsSplit <- split(allHapBlocks, allHapBlocks$Chemical)

haploidPlotter <- lapply(allHapsSplit, function(chem) {
    fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/", chem$Chemical[1], "_50snp_window_unfiltered")
    chemDFs <- split(chem, chem$id)
    hapPlot <- plotting.haploid.haplotypes(chemDFs, fileName)
} )

##################################  Plotting RHC haplotype blocks single chromosome  #####################
newFounderNames <- read.table("newFounderNames.txt", sep = '\t', colClasses = 'character')
newFounderNames <- newFounderNames[,1]
#newFounderNames[9] <- 'A12'
#newFounderNames[11] <- 'B6'
#newFounderNames[15] <- 'B12'

bas01RhcPlotChr2 <- plotting.zoomed.in.haploid.haplotypes(definedHapBlocks[c(1:9)], '18F13v1 recombinant haploid clones chrX', 10)
bas02RhcPlotChr2 <- plotting.zoomed.in.haploid.haplotypes(definedHapBlocks[10:19], '18F13v2 recombinant haploid clones chrX', 10)



# ============================================================================
# Trouble-shooting

rhc_interval_dfs <- splitDFs
fname <- fileName

plotting.haploid.haplotypes <- function(rhc_interval_dfs, fname) {
    offsets <- read.table("newoffsets.txt", header = T)
    offsets$lines <- ceiling(offsets[,2]/50)
    offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
    offsets$chr <- 1:17
    offsets <- offsets[-nrow(offsets),] # no mitochondira
    g_l <- c(0, cumsum(offsets$len))
    ch.bounds <- c(0, g_l[1:16] + offsets[,2])
    #correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
    max.pos <- g_l[length(g_l)]
    mid.ch <- diff(ch.bounds)/2
    midpt.ch <- ch.bounds[2:17] - mid.ch
    plotting_factors <- c(topHaps, 'unknown', 'unresolvable')
    plotting_colors <- c(mycols, 'darkgray')
    color_df <- data.frame(factors = plotting_factors, colors = plotting_colors, stringsAsFactors = FALSE)
    title <- strsplit(fname, "/")[[1]]
    pTitle <- title[[length(title)]]
    changeRle <- rle(filteredDF$Chemical)
    changeRle2 <- cumsum(changeRle$lengths)
    changeIdces <- changeRle2+ 1
    changeIdces <- changeIdces[-length(changeIdces)]
    changeIds <- filteredDF$id[changeIdces]
    pdf(paste0(fname, '.pdf'), height = 4, width = 6)
    par(mar = c(0,0,0,6), oma = c(3,2,3,2), xpd = TRUE)
    par(mgp = c(1.7,0.75,0))
    layout(matrix(1:(length(rhc_interval_dfs)+3), length(rhc_interval_dfs)+3, byrow = TRUE))
    counts <- 0
    plot_all <- lapply(rhc_interval_dfs, function(df) {
        counts <<- counts + 1
        if(df$Chemical[1] == "base") {
            id <- "base"} else {
                id <- paste0("wk11\n", df$id[1])
            }
        #print(id)
        #flush.console()
        df$mfh[df$mfh == '?'] <- 'unresolvable' 
        #df$mfhp <- df$mfh
        df$mfh <- gsub(",", "", df$mfh)
        df$mfh[df$mfh != 'unresolvable' & !df$mfh %in% color_df$factors] <- 'unknown'
        df$cols <- unlist(lapply(df$mfh, function(x) {
            color_df$colors[color_df$factors == x]
        } ) )
        #par(mar = c(1.25,0,0,0), mgp = c(1.2,1,0), oma = c(2,4,3,1))
        #par(mar = c(0,0,0,0), oma = c(1,0.25,0,0.5), mgp = c(0.75,0.5,0))
        xmin <- 0
        xmax <- max.pos
        if(df$id[1] %in% changeIds) {
            plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(xmin, xmax), ylim= c(0, 1), axes="FALSE", yaxs="i", xaxs="i")
            #title(paste0("week 11 ", df$Chemical[1]), cex.main = 0.5, line = 0, xpd = NA)
        } 
        plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(xmin, xmax), ylim= c(0, 1), axes="FALSE", yaxs="i", xaxs="i")
        plotting <- lapply(seq(1,nrow(df), 2), function(x){
            rect(df$gp[x], 0.1, df$gp[x+1], 0.9, density = NA, col = df$cols[x], border = NA)
        } )
        #axis(2, at = 0.5, labels = id, las = 1, cex.axis = 0.5, lwd = 0, line = -0.5)
        clip(0,xmax,0,1)
        #abline(h= 0, col = 'black', lwd = 1)
        abline(v= ch.bounds, col = 'black', lwd = 0.5, lty = 1)
        if(counts == 1) {
            abline(h = 1, col = 'black', lwd = 1)
            #title(main = "base", cex.main = 0.5, line = 0.25, xpd = NA)
            legend(grconvertX(1.005, "npc", "user"),grconvertY(1.25, "npc", "user"), legend = unique(plotting_factors), col = plotting_colors, lwd = 2, cex = 0.75, xpd = NA, bty = "n") }
        if(counts == length(rhc_interval_dfs)) {
            axis(1, at = midpt.ch, labels = as.roman(1:16), las = 1, cex.axis = 0.75, lwd = 0, line = -0.5)
            clip(0,xmax,0,1)
            abline(h= 0, col = 'black', lwd = 1)}
    } )		
    dev.off()
}


allHapsWant <- allHapBlocks[allHapBlocks$week %in% c(0, 11),]
fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/FigureS_")
chemDFs <- split(allHapsWant, allHapsWant$id)
filteredDFs <- lapply(chemDFs, function(x) {
    unresolvLength <- sizing.unknown.blocks(x)
    fracUnresolv <- unresolvLength/max.pos
    if(fracUnresolv < 0.25) {return(x)}
})
filteredDF <- do.call(rbind, filteredDFs)
splitDFs <- split(filteredDF, filteredDF$id)
hapPlot <- plotting.haploid.haplotypes(splitDFs, fileName)


haploidPlotter <- lapply(allHapsSplit, function(chem) {
    fileName <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Haploid_clones/", chem$Chemical[1], "_50snp_window_filtered_notitle")
    chemDFs <- split(chem, chem$id)
    filteredDFs <- lapply(chemDFs, function(x) {
        unresolvLength <- sizing.unknown.blocks(x)
        fracUnresolv <- unresolvLength/max.pos
        if(fracUnresolv < 0.25) {return(x)}
    })
    filteredDF <- do.call(rbind, filteredDFs)
    splitDFs <- split(filteredDF, filteredDF$id)
    hapPlot <- plotting.haploid.haplotypes(splitDFs, fileName)
} )
