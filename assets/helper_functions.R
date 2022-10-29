##############################################################################

# Title:       Helper functions
# Author:      Robert Linder
# Date:        2022
# Description: For providing functions used in multiple parts of this analysis pipeline.
 
##############################################################################

# ============================================================================
# Function details     

# position.offsetrer:  A detailed description

# ============================================================================
# Custom functions

add.alpha <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
              rgb(x[1], x[2], x[3], alpha=alpha))  
}

calcstats_asin.sqrt.hap.avg.dif.transforms <- function(treatmentDT, treatCol, cutoffConserv, cutoffLib) {
    tictoc::tic('total')
    id <- treatmentDT[, get(treatCol)]
    treatmentDT[, id := gsub("R[0-9][0-9]$", "", id)]
    cat(paste0(treatmentDT$id[1], '\n'))
    flush.console()
    reps <- treatmentDT[, uniqueN(Replicate)]
    treatDifs <- treatmentDT[, tstrsplit(asinSqrtFreqDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
    founderDT <- treatmentDT[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
    founderFreqs <- treatmentDT[, tstrsplit(baseRawFreqsCol, split = ";", type.convert = TRUE, fixed = TRUE)]
    founderDT$gp <- treatmentDT$gp
    treatDifs$gp <- treatmentDT$gp
    founderFreqs$gp <- treatmentDT$gp
    avgDT <- treatDifs[, .(avgDifs = unlist(lapply(.SD, mean)), varDifs = unlist(lapply(.SD, var))), .SDcols = c(paste0("V", 1:17)), by = gp]
    foundersDT <- founderDT[, .(founders = unlist(lapply(.SD, function(x) {x[1]}))), .SDcols = c(paste0("V", 1:17)), by = gp]
    avgDTbase <- founderFreqs[, .(baseFreqs = unlist(lapply(.SD, function(x) {x[1]}))), .SDcols = c(paste0("V", 1:17)), by = gp]
    avgAllDT <- cbind(avgDT, avgDTbase$baseFreqs)
    kCol <- avgAllDT[V2 > cutoffConserv, .N, by = gp]
    statDT <- avgAllDT[V2 > cutoffConserv, .(sumHapDifs = sum(avgDifs^2), sumVarDifs = sum(varDifs)), by = gp][, k := kCol[, N]]
    sigma_a <- statDT[, sum(sumVarDifs)/sum(k)]
    sigma_b <- (0.004)^2 + (0.01)^2
    V <- sigma_a/reps + sigma_b
    kColLib <- avgAllDT[V2 > cutoffLib, .N, by = gp]
    statLibDT <- avgAllDT[V2 > cutoffLib, .(sumHapDifs = sum(avgDifs^2), sumVarDifs = sum(varDifs)), by = gp][, k := kCol[, N]]
    lodDT <- statLibDT[, C := sumHapDifs/V][, LOD := -pchisq(C,k,lower.tail=FALSE,log.p=TRUE)/log(10)]
    infoDT <- unique(treatmentDT, by = "gp")
    lodAllDT <- lodDT[infoDT, on = .(gp), nomatch = 0][, c("gp", "chr", "pos"
                                                           , "Chemical", "id", "Week", "chemWeek", "sumHapDifs", "sumVarDifs", "k", "C", "LOD"), with = FALSE
    ][, c("sigma_a", "V") := .(sigma_a, V)]
    tictoc::toc()
    lodAllDT
}

calcstats_ind.hap.dif.heterozygosities <- function(treatmentDT, baseDT, treatCol) {
    tictoc::tic('total')
    id <- treatmentDT[, get(treatCol)]
    cat(paste0(id[1], '\n'))
    flush.console()
    baseHeterozygosityCol <- baseDT$heterozygosity
    #hetSqDevsCol <- sum((sort(baseDT$heterozygosity)-sort(treatmentDT$heterozygosity))^2)/nrow(baseDT)
    hetAvgDevCol <- sum(baseDT$heterozygosity - treatmentDT$heterozygosity)/nrow(baseDT)
    hetUpCol <- length(which(baseDT$heterozygosity - treatmentDT$heterozygosity < 0))/nrow(baseDT)
    baseIDCol <- baseDT$id
    fullDT <- cbind(treatmentDT, baseHeterozygosityCol, hetAvgDevCol, hetUpCol, baseIDCol)
    tictoc::toc()
    fullDT
}

calcstats_ind.hap.dif.transforms <- function(treatmentDT, baseDT, treatCol) {
    tictoc::tic('total')
    id <- treatmentDT[, get(treatCol)]
    cat(paste0(id[1], '\n'))
    flush.console()
    baseDT <- baseDT[grepl(treatmentDT$chemWeek[1], baseDT$id)]
    treatFreqs <- treatmentDT[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
    asinSqrtTreatFreqs <- asin(sqrt(treatFreqs))
    baseFreqs <- baseDT[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
    asinSqrtBaseFreqs <- asin(sqrt(baseFreqs))
    freqDifs <- treatFreqs - baseFreqs
    difSign <- freqDifs[, lapply(.SD, function(x) gsub("^-.*", "Decreasing", x))
    ][, lapply(.SD, function(x) gsub("^0$", "Neutral", x))
    ][, lapply(.SD, function(x) gsub("^([1-9].*|0.0*\\d+|0.\\d+)", "Increasing", x))]
    asinSqrtFreqDifs <- asinSqrtTreatFreqs - asinSqrtBaseFreqs
    freqDifsSq <- freqDifs^2
    asinSqrtTreatFreqsCol <- unite(asinSqrtTreatFreqs, col = "asinSqrtTreatFreqs", sep = ";", remove = TRUE, na.rm = TRUE)
    asinSqrtBaseFreqsCol <- unite(asinSqrtBaseFreqs, col = "asinSqrtBaseFreqsCol", sep = ";", remove = TRUE, na.rm = TRUE)
    freqDifsCol <- unite(freqDifs, col = "freqDifs", sep = ";", remove = TRUE, na.rm = TRUE)
    freqDifsSqCol <- unite(freqDifsSq, col = "freqDifsSq", sep = ";", remove = TRUE, na.rm = TRUE)
    sumFreqDifsSqCol <- rowSums(freqDifsSq, na.rm = TRUE)
    difSignCol <- unite(difSign, col = "difSign", sep = ";", remove = TRUE, na.rm = TRUE)
    asinSqrtFreqDifsCol <- unite(asinSqrtFreqDifs, col = "asinSqrtFreqDifs", sep = ";", remove = TRUE, na.rm = TRUE)
    baseRawFreqsCol <- baseDT$collapsedFreqs
    baseHeterozygosityCol <- baseDT$heterozygosity
    hetSqDevsCol <- sum((sort(baseDT$heterozygosity)-sort(treatmentDT$heterozygosity))^2)/nrow(baseDT)
    baseIDCol <- baseDT$id
    fullDT <- cbind(treatmentDT, freqDifsCol,freqDifsSqCol, sumFreqDifsSqCol, difSignCol, asinSqrtTreatFreqsCol, asinSqrtBaseFreqsCol, asinSqrtFreqDifsCol, baseRawFreqsCol, baseHeterozygosityCol, hetSqDevsCol, baseIDCol)
    tictoc::toc()
    fullDT
}

calcstats_snp.dif.transforms <- function(treatmentSample, baseSample, mathFUN, snpDT) {
    treatmentDT <- snpDT[, names(snpDT) == treatmentSample, with = FALSE]
    baseDT <- snpDT[, names(snpDT) == baseSample, with = FALSE]
    mathFun <- get(mathFUN)
    difDT <- mathFun(treatmentDT - baseDT)
    calcStatDT <- data.table()
    calcStatDT[, c("chr", "POS", "gp", "id", paste0(treatmentSample, '-', baseSample)) := .(snpDT$chr, snpDT$POS, snpDT$gp, paste0(treatmentSample, '-', baseSample), difDT[[1]])]
    calcStatDT
}

covcalc_create.chr.cov.dts <- function(window, samplePattern) {
    function(snpCovDF) {
        tic('total')
        covListDTs <- lapply(unique(names(snpCovDF)[grepl(samplePattern, names(snpCovDF))]), function(f) {
            print(f)
            tic('this cycle')
            #print(names(snpCovDT))
            #flush.console()
            sample <- snpCovDF[names(snpCovDF) == f]
            data <- cbind(snpCovDF[1:3], sample)
            names(data) <- c(names(data)[1:3], f)
            data <- data[data$chr != 17,]
            rdna <- which(data[,1] ==12 & data[,2] >= 451418 & data[,2] <= 468931)
            if(length(rdna) > 0) {data <- data[-rdna,]} else {
                data <- data}
            split.data <- split(data, data[,1])
            win.measures <- lapply(split.data, function(i) {
                wins <- seq(min(i[,2]), max(i[,2]), window)
                m <- sapply(wins, function(j) {
                    interval <- j:(j+window-1)
                    positions <- which(i[,2] %in% interval)
                    if(length(positions) > 0) {
                        avg <- mean(i[,4][i[,2] %in% interval], na.rm = T) } else {
                            avg <- NA} 
                    return(avg) 
                } )
                cbind(c = i[,1][1], w = as.vector(wins), m = as.vector(m)) 
            } )
            data <- data.frame(do.call("rbind", win.measures), stringsAsFactors = F)
            data$cp <- apply(data[,1:2], 1, paste, collapse = "_")
            data <- data[!(data$cp %in% repeats$cp),]
            mean_cov <- mean(data[,3], na.rm = TRUE)
            data$relCov <- data[,3]/mean_cov
            data$gp <- data$w + g_l[data$c]
            #data <- data[!(data$gp %in% badSnpsDT$gp),]
            names(data)[3] <- 'avgCov'
            names(data)[2] <- 'POS'
            names(data)[1] <- 'chr'
            data$id <- rep(f, nrow(data))
            data$gwideCov <- rep(mean_cov, nrow(data))
            data <- data[, c(1,2,6,4,7,3,5,8)] 	
            toc()
            data <- as.data.table(data)
            data
        } )
    }
}

covcalc_create.gwide.cov.dts <- function(snpCovDT, samplePattern) {
    tic('total')
    covListDTs <- lapply(unique(names(snpCovDT)[grepl(samplePattern, names(snpCovDT))]), function(f) {
        print(f)
        tic('this cycle')
        #print(names(snpCovDT))
        flush.console()
        setDTthreads(1)
        subsetDT <- snpCovDT[, c('chr', 'POS', 'gp', f), with = FALSE 
        ]
        subsetDT2 <- subsetDT[!which(is.na(subsetDT[, f, with = FALSE]))
        ][!.(12, 451418:468931), on = c('chr', 'POS')
        ][!.(17), on = 'chr'
        ][, c('cp', 'gp', 'id') := .(paste(chr, POS, sep = "_"), POS + g_l[chr], rep(f, .N))
        ][!cp %in% repeats$cp
        ][order(gp)
        ][, .(id = id[1], avgCov = mean(get(f), na.rm = T))]
        toc()
        subsetDT2
    } )
    toc()
    do.call(rbind, covListDTs)
}

covcalc_create.norm.cov.dts <- function(treatmentDT, baseDT) {
    tictoc::tic('total')
    cat(paste0(treatmentDT$id[1], '\n'))
    flush.console()
    joinedDTs <- treatmentDT[baseDT, on = .(gp), nomatch = 0]
    calcCovDifsDT <- joinedDTs[, normCov := (avgCov/i.avgCov)*(i.gwideCov/gwideCov)]
    tictoc::toc()
    calcCovDifsDT
}
covcalc_filter.low.cov.samples <- function(covCutoff) {
    function(gwideCovDT) {
        gwideCovDT[avgCov < covCutoff, .(id)]
    }
}

filesplitter_split.collapsed.haps.basic.samples <- function(listOfTreatments, idCol) {
    function(hapDT) {
        indTxDTs <- lapply(listOfTreatments, function(x) {
            tic()
            print(x)
            flush.console()
            txDT <- hapDT[get(idCol) == x
            ][, 'id' := NULL
            ][, c('cp', 'id') := .(paste(chr, pos, sep = "_"), get(idCol))]
            if(grepl("CNTL", txDT$id[1])) {
                chem <- txDT$id[1]} else {
                    chem <- paste(strsplit(txDT$id[1], "")[[1]][4:10], collapse = "")
                    if(grepl("NA", chem)) {chem <- txDT$id[1]} else {chem <- chem}
                }
            if(grepl("BAS02", chem)) {
                txDT$Chemical <- treatKeyDT[V1 == chem, Chemical]
                txDT$Week <- 0
            } else {
                txDT$Chemical <- treatKeyDT[V1 %like% chem, Chemical][1]
                txDT$Week <- paste(strsplit(txDT$id[1], "")[[1]][4:5], collapse = "")
            }
            txDT$Replicate <- paste(strsplit(txDT$id[1], "")[[1]][14:16], collapse = "")
            txDT$chemWeek <- paste(txDT$Chemical, txDT$Week, sep = "_")
            if(grepl("H[0-9][0-9]$", x)) {
                txDT$haploid <- substr(x, nchar(x)-1, nchar(x))
            }
            txDT2 <- txDT[, c(idCol) := NULL][chr != 17]
            if("founderfreqs" %in% names(txDT2)) {
                sepCuts <- lapply(txDT2$cutree, function(x) as.numeric(strsplit(x, ";", fixed = TRUE)[[1]]))
                sepFreqs <- lapply(txDT2$founderfreqs, function(x) as.numeric(strsplit(x, ";", fixed = TRUE)[[1]]))
                pooledFreqs <- lapply(1:length(sepCuts), function(x) {
                    mergeFreqs <- tapply(sepFreqs[[x]], sepCuts[[x]], sum)
                    sumFreqs <- sum(mergeFreqs)
                    het <- 1 - sum(mergeFreqs^2)
                    mergedFreqs <- paste(mergeFreqs, collapse = ";")
                    mergeFounders <- paste(tapply(founderNames, sepCuts[[x]], paste, collapse = ""), collapse = ";")
                    return(c(mergedFreqs, mergeFounders, sumFreqs, het))
                } )
                newFounderCols <- as.data.table(do.call(rbind, pooledFreqs))
                names(newFounderCols) <- c("collapsedFreqs", "collapsedFounders", "sumFreqs", "heterozygosity")
                txDT2 <- cbind(txDT2, newFounderCols)
                txDT2 <- txDT2[, 'founderfreqs' := NULL]
                txDT2$heterozygosity <- as.numeric(txDT2$heterozygosity)
                txDT2$sumFreqs <- as.numeric(txDT2$sumFreqs)
            }
            setcolorder(txDT2, c('chr', 'pos', 'gp', 'cp', 'Chemical', 'Replicate', 'id'))
            toc()
            txDT2
        } )
    }
}

hapfilereform_check.haplotype.file <- function(formattedHapDT, idCol) {
    sampleCountsDT <- formattedHapDT[, .(length = .N), by = .(sample = get(idCol))]
    if(length(rle(sampleCountsDT[, length])$values) != 1) {stop("Not all samples have the same number of haplotype calls.")} else {
        print("Good to go!")
        flush.console()}
}

hapfilereform_reformat.haplotype.file <- function(samplePattern, basePop, colsKeeping, idCol) {
    function(hapDT) {
        hapDT2 <- hapDT[get(idCol) %like% samplePattern | get(idCol) %like% basePop, .SD, .SDcols = colsKeeping
        ][,chr := as.integer(as.roman(substr(chr, 4, nchar(chr))))
        ][chr == 1000, chr := 17
        ][order(chr, pos, get(idCol))
        ][,gp := pos + g_l[chr]
        ][order(get(idCol), chr, pos)
        ][get(idCol) %like% "_", c(idCol) := gsub("_", "", get(idCol))
        ][, id := paste(get(idCol), chr, pos, gp, sep = "_")
        ]
        hapDT2
    }
}

posoff_chr.bounds <- function(offsets_dt) {
	c(0, cumsum(offsets_dt$len))
 }

reading_read.in.dts.as.list <- function(projectRootDir, folderPath, analysisType, samplePattern) {
    files <- list.files(path = paste0(projectRootDir, folderPath), pattern = paste0(analysisType, "_DT.txt$"))
    files <- files[grepl(samplePattern, files)]
    print(files)
    flush.console()
    lapply(files, function(x) {
        filePath <- paste0(projectRootDir, folderPath, x)
        fread(filePath, header = T)
    } )
}

rep.list <- function(object, repObject) {
    rep(list(object), length(repObject))
}

test.match.order <- function(x,y) {
    if (isTRUE(all.equal(x,y))) print('Perfect match in same order')
    if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in wrong order')
    if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) print('No match')
}

writing_write.file <- function(projectRootDir, folderPath, fileName) {
    function(DT) {
        print(fileName)
        flush.console()
        fileName <- paste0(fileName, "_DT.txt")
        filePath <- paste0(projectRootDir, folderPath, fileName)
        fwrite(DT, filePath, sep = "\t")
    }
}

writing_write.lists.of.dts.to.txt.files <- function(projectRootDir, folderPath, analysisType, idCol) {
    function(listOfDTs) {
        lapply(listOfDTs, function(x) {
            fileName <- paste0(x[1, get(idCol)], "_", analysisType, "_DT.txt") 
            filePath <- paste0(projectRootDir, folderPath, fileName )
            fwrite(x, filePath, sep = "\t")
        } )
    }
}
 

 
# ============================================================================
# Trouble-shooting
