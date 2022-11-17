#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Finding aneuploidies
# Description: Finds aneuploidies and infers their bounds using a Hidden Markov Model by using normalized coverage relative to the base population. 

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("Usage: analysis_step14-Aneuploidies_calling_HMM.R <normalized_coverages> <low_cov> <offsets> <helper>", call.=FALSE)
}

norm_covs <-  args[1]
low_cov <- args[2]
offsets <- args[3]
helper <- args[4]

# ============================================================================
# Load packages and sourced files

library(data.table)
library(tictoc)
library(HMM)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()
options(digits = 10)

# ============================================================================
# Custom functions

flipstate <- function(t) {
    t <- split(t, t$chr)
    t <- lapply(t, function(d) {
        nas <- which(is.na(d$h))
        #print(d[,1][1])
        newvals <- sapply(nas, function(x) {
            #print(x)
            up <- na.omit(rev(d$h[1:(x-1)]))[1]  ### looking at call next to (to the left of) the NA in $h
            down <- na.omit(d$h[(x+1):length(d$h)])[1] ### looking at call to the right of the NA in $h
            if (x != 1 & x != length(d$h)) {
                if (is.na(up) == FALSE & is.na(down) == FALSE) {	### checking that no NAs to left or right of NA currently looking at
                    if (up == down) {return(up)}	### if have same call to left and right of NA, turn NA into that call
                    else {return(NA)}					
                } else if (is.na(up) == FALSE & is.na(down) == TRUE) {	### use call to the left of NA
                    return(up)
                } else if (is.na(up) == TRUE & is.na(down) == FALSE) {	### use call to the right of NA
                    return(down)
                }
            } else if (x == 1) {	### dealing with an NA as the first call
                return(na.omit(d$h[(x+1):length(d$h)])[1])	### use call directly to the right
            } else if (x == length(d$h)) {	### dealing with an NA as the last call
                return(na.omit(rev(d$h[1:(x-1)]))[1])	### use call directly to the left
            }
            
        })
        if (length(nas) > 0) {d$h[nas] <- unlist(newvals)} ### replacing all NAs in $h with the corrected values (0 or 1) called above
        return(d)
    })
    t <- do.call("rbind", t)
    return(t)
}

callgeno <- function(t) {
    fname <- t$treatment[1]
    print(fname)
    flush.console()
    t <- t[t$chr %in% 1:16,]
    na.data <- t[is.na(t$Sign),]
    #u <- round(mean(t[,3]), digits = 2)
    t <- na.omit(t)
    #cond <- (t$f >= 0.75 | t$f <= 0.25)
    #bt.data <- t[cond == FALSE,]
    #t <- t[cond,]
    # new hmm params:  initHMM(as.character(0:1), as.character(0:1), c(.5,.5), matrix(c(.8,.01,.2,.99),2),matrix(c(.6,.4,.4,.6),2))
    hmm <-initHMM(as.character(0:1), as.character(0:1), c(.7,.3), matrix(c(.8,.01,.2,.99),2),matrix(c(.6,.4,.4,.6),2)) ### normally use
    t2 <- split(t, t$chr)
    t2 <- lapply(t2, function(d) {
        d$h <- unlist(as.numeric(viterbi(hmm, as.character(d$Sign))))
        d
    })
    t2 <- do.call("rbind", t2)
    if (nrow(na.data) > 0) {
        na.data$h <- NA
        t2 <- rbind(t2,na.data)
    }
    #if (nrow(bt.data) > 0) {
    #bt.data$h <- NA
    #t <- rbind(t,bt.data)
    #}
    t2 <- t2[order(t$gp),]
    t2 <- flipstate(t2)
    t3 <- na.omit(t2)
    #cond1 <- (t2[,9] >= 0.75 | t2[,9] <= 0.25) 
    #t2 <- t2[cond1,]
    png(paste0(filePath,fname,"_updated3.png"),  width = 10, height = 5, units = 'in', res = 150)
    par(mar = c(2,0,1,0), oma = c(1,7,3,4), mgp = c(1.2,1,0))
    layout(matrix(1:16,1), widths = c((ch.sizes[-17])/max.pos[16]))
    plotter <- sapply(1:16,function(ch) {
        plot(t3$POS[t3$chr == ch], t3$normCov[t3$chr == ch], col = mypal[ch] , pch = 20, xaxt = "n", xaxs = "i", bty = "n", xlab = "", ylim = c(0,3), yaxt = "n")
        points(t3$POS[t3$chr == ch], t3$h[t3$chr == ch] + 1, type = "l", col = mypal[ch], lwd = 1.5) 
        if(ch == 1) {
            mtext("Normalized coverage", side = 2, line = 4)
            axis(2, at = seq(0,3,0.5), labels = seq(0,3,0.5), cex.axis = 1.5, las = 1)}
        if(ch == 10) {
            title(main = fname, xpd = NA)} 
        if(ch == 16) {mtext(side = 4, line = 2, paste("Average coverage =", round(t3$gwideCov[1], 2)))}
        axis(1, at = ch.sizes[ch]/2, labels = as.roman(ch), cex.axis = 1.5, tick = T) } )
    dev.off()	
    return(as.data.table(t2))
}

find.significant.regions <- function(df) {
    tic('total')
    print(df$id[1])
    flush.console()
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    windows <- df[df$h == 1,]
    sig_rows <- as.numeric(rownames(windows))
    if(nrow(windows) > 0) {
        if(nrow(windows) > 1) {
            return_boundaries <- unlist(lapply(1:length(sig_rows[-length(sig_rows)]), function(x) {
                start <- diff(c(sig_rows[x], sig_rows[x + 1]))
                if(start > 1 & x == 1 & length(sig_rows) == 2) {return(c(rep(sig_rows[x], 2), rep(sig_rows[x+1], 2)))
                } else if(start > 1 & x == 1) {return(c(rep(sig_rows[x], 2), sig_rows[x+1]))
                } else if(start == 1 & x == 1 & x + 1 == length(sig_rows)) {return(c(sig_rows[x], sig_rows[x+1]))
                } else if(start == 1 & x == 1) {return(sig_rows[1])
                } else if(start > 1 & x == length(sig_rows) - 1) {return(c(sig_rows[x], rep(sig_rows[x+1], 2)))
                } else if(start == 1 & x + 1 == length(sig_rows)) {return(sig_rows[length(sig_rows)])
                } else if(start > 1) {return(c(sig_rows[x], sig_rows[x+1]))
                } else if(start == 1 & windows$chr[x] != windows$chr[x + 1]) {return(c(sig_rows[x], sig_rows[x+1]))
                }
            } ) ) 
        } else {return_boundaries <- c(1,1)}
        windows_df <- df[return_boundaries,]
        rownames(windows_df) <- 1:nrow(windows_df)
        if(nrow(windows_df) > 2) {
            too_close <- unlist(lapply(seq(2, nrow(windows_df) - 2, 2), function(x) {
                if(windows_df$chr[x+1] == windows_df$chr[x]) {
                    distance <- windows_df$gp[x + 1] - windows_df$gp[x]
                    if(distance < 10000) {return(x)} ### originally 50kb
                } } ) )
            collapsing_windows <- c(too_close, too_close + 1)
            if(length(collapsing_windows) > 1){
                windows_df <- windows_df[-c(collapsing_windows),]} else{
                    windows_df <- windows_df} }
        rownames(windows_df) <- 1:nrow(windows_df)
        if(nrow(windows_df) %% 2 != 0) {break}
        toc()
        windows_df
    } 
}

find.most.significant.subregions <- function(windows, df) {
    windows <- as.data.frame(windows, stringsAsFactors = FALSE)
    windows$avgNormCov <- 0
    windows$maxStartPos <- 0
    windows$maxEndPos <- 0
    windows$maxStartGP <- 0
    windows$maxEndGP <- 0
    windows$maxValue <- 0
    sigCols <- "normCov"
    windows1.2 <- lapply(seq(1,nrow(windows), 2), function(x) {
        print(x)
        flush.console()
        peakpos <-  df[df$gp >= windows$gp[x] & windows$gp[x+1] >= df$gp,]
        meanPeak <- mean(peakpos[,sigCols], na.rm = TRUE)
        if(meanPeak > 1) {
            maxPeak <- max(peakpos[,sigCols], na.rm = TRUE)
            sig_rows <- which(peakpos[, sigCols] == maxPeak)
            rows <- c(min(sig_rows), max(sig_rows))} else {
                minPeak <- min(peakpos[,sigCols], na.rm = TRUE)
                sig_rows <- which(peakpos[, sigCols] == minPeak)
                rows <- c(min(sig_rows), max(sig_rows))}
        windows$avgNormCov[x:(x+1)] <- meanPeak
        windows$maxValue[x:(x+1)] <- peakpos$normCov[rows[1]]
        windows$maxStartPos[x:(x+1)] <- peakpos$POS[rows[1]]
        windows$maxStartGP[x:(x+1)] <- peakpos$gp[rows[1]]
        if(rows[1] == rows[2]) {
            windows$maxEndPos[x:(x+1)] <- peakpos$POS[rows[2]] + 0.5
            windows$maxEndGP[x:(x+1)] <- peakpos$gp[rows[2]] + 0.5} else {
                windows$maxEndPos[x:(x+1)] <- peakpos$POS[rows[2]]
                windows$maxEndGP[x:(x+1)] <- peakpos$gp[rows[2]]}
        windows[x:(x+1),]
    } )
    windows1.3 <- do.call(rbind, windows1.2)
    windows1.3 <- as.data.table(windows1.3)
    windows1.3[avgNormCov < 1, Type := "Deletion"]
    windows1.3[avgNormCov > 1, Type := "Duplication"]
    windows1.3
}

test.match.order <- function(x,y) {
    if (isTRUE(all.equal(x,y))) print('Perfect match in same order')
    if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in wrong order')
    if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) print('No match')
}

# ============================================================================
# Load data

NormCovDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "/Treatment_vs_base_coverage_tables_2kb/", analysisType = "normalized_coverage", samplePattern = ".*")
lowCovDTs <- fread(low_cov, header = T)
offsets <- read.table(offsets, header = T)
offsets$lines <- ceiling(offsets[,2]/50)
offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
offsets$chr <- 1:17
g_l <- c(0, cumsum(offsets$len))
ch.bounds <- c(0, g_l[1:17] + offsets[,2])
ch.sizes <- diff(ch.bounds)
#correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
max.pos <- cumsum(ch.sizes)
mid.ch <- diff(ch.bounds)/2
midpt.ch <- ch.bounds[2:18] - mid.ch

filePath <- paste0(projectDir, "/Plots/Ind_normalized_coverage_mutations_plots/")
mypal <- rep(RColorBrewer::brewer.pal(8,'Dark2'),2)

# ============================================================================
# Run an HMM to find aneuploidies and large duplications/deletions

allCovDTs <- do.call(rbind, NormCovDTs)
allCovDTs$treatment <- unlist(lapply(1:nrow(allCovDTs), function(x) {
    split <- paste(strsplit(allCovDTs$id[x], "")[[1]][9:16], collapse = "")
    #print(split)
    #flush.console()
    if(split == "YP") {return("control")} else{
        return(split)}
} ) )
stillGoodDT <- allCovDTs[!id %in% lowCovDTs$id][!id %like% "000YP"][!id %like% "CNTL"]
# old params 1.3, 0.7
stillGoodDT[, Sign := as.numeric(normCov >= 1.25 | normCov <= 0.75)]
stillGoodDF <- as.data.frame(stillGoodDT, stringsAsFactors = FALSE)
splitDF <- split(stillGoodDF, stillGoodDF$id)

runHMM <- Map(callgeno, splitDF)

savingNormalizedHMMCovsDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "/Tables/Treatment_vs_base_coverage_tables_2kb/", analysisType = "normalized_coverage_hmm", idCol = "id")
writeNormalizedCovsFiles <- savingNormalizedHMMCovsDTs(runHMM)

# ============================================================================
#  Putting together windows of aneuploidies

aneuploidWins <- Map(find.significant.regions, runHMM)
aneuploidWinsDTs <- Filter(Negate(function(i) is.null(unlist(i))), aneuploidWins)

orderCovs <- unlist(lapply(aneuploidWinsDTs, function(x) return(x$id[1])))
orderSplitDF <- unlist(lapply(splitDF, function(x) return(x$id[[1]])))
notIn18 <- which(!orderCovs %in% orderSplitDF)
notInnext18 <- which(!orderSplitDF %in% orderCovs)
if(length(notIn18) > 0) {aneuploidWinsDTs <- aneuploidWinsDTs[-notIn18]} else{aneuploidWinsDTs <- aneuploidWinsDTs}
if(length(notInnext18) > 0) {splitDF <- splitDF[-notInnext18]} else{splitDF <- splitDF}
newOrderCovs18 <- unlist(lapply(aneuploidWinsDTs, function(x) return(x$id[1])))
newOrderCNVs18<- as.character(unlist(lapply(splitDF, function(x) return(x$id[[1]]))))
checkOrder18 <- test.match.order(newOrderCovs18, newOrderCNVs18)

aneuploidWinsComplete <- Map(find.most.significant.subregions, aneuploidWinsDTs, splitDF)
aneuploidWinsCompleteDT <- Filter(Negate(function(i) is.null(unlist(i))), aneuploidWinsComplete)

savingNormalizedHMMCovsDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "/Tables/Aneuploid_windows/", analysisType = "aneuploid_windows", idCol = "id")
writeNormalizedCovsFiles <- savingNormalizedHMMCovsDTs(aneuploidWinsCompleteDT)

# ============================================================================
# Trouble-shooting