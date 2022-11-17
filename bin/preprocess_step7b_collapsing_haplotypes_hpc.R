#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Collapsing haplotypes
# Description: Collapses founder haplotypes that are genetically identical at certain genomic locations; this is done on a per position basis, creating synthetic founder haplotypes that are a fusion of the founder strains indistinguishable at a particular genomic locatin. 

##############################################################################
# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    stop("Usage: preprocess_step7b_collapsing_haplotypes_hpc.R <haplotypes> <chemical> <treatment_key> <founders> <helper>", call.=FALSE)
}

haplotypes <- args[1]
chem <-  args[2]
treatment_key <- args[3]
founders <- args[4]
helper <- args[5]

# ============================================================================
# Load packages and sourced files

library(data.table)
library(tictoc)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

merge.master <- function(M) {
    counter <<- counter + 1
    setTxtProgressBar(pb, counter)
    test = M[[1]]
    for(i in 2:length(M)){
        newtest1 = merge.down(test, M[[i]])
        change="Yes"
        while(change=="Yes"){
            newtest2 = simplify(newtest1)
            newtest3 = merge.across(newtest2)
            #cat(i,"\t",newtest2,"\t",newtest3,"\n")
            if(length(newtest2) == length(newtest3) & sum(newtest2 == newtest3) == length(newtest2)){
                change="No"
                test=newtest3
            }else{
                newtest1 = newtest3
            }
        }
    }
    test
}

merge.down=function(A,B){
    out=A
    for(i in 1:length(A)){
        for(j in 1:length(B)){
            if(sum(unlist(strsplit(A[i],'')) %in% unlist(strsplit(B[j],''))) + sum(unlist(strsplit(B[j],'')) %in% unlist(strsplit(A[i],''))) != 0){
                out[i] = paste(out[i],B[j],sep='')
            }
        }
    }
    out
}

simplify=function(A){
    out = rep(NA,length(A))
    for(i in 1:length(A)){
        out[i] = paste(unique(sort(unlist(strsplit(A[i],'')))),collapse='')
    }
    out
}

merge.across=function(A){
    L = length(A)
    out = A
    i = 1
    while(i <= L-1){
        j = i+1
        keep=rep(1,length(out))
        while(j <= L){
            if(sum(unlist(strsplit(out[i],'')) %in% unlist(strsplit(out[j],''))) + sum(unlist(strsplit(out[j],'')) %in% unlist(strsplit(out[i],''))) != 0){
                out[i] = paste(out[i],out[j],sep='')
                keep[j] = 0
            }
            j=j+1
        }
        out=out[as.logical(keep)]
        i=i+1
    }
    out
}

# ============================================================================
# Load data

hapFreqs <- fread(haplotypes, header = T)
treatKeyDT <- fread(treatment_key, header = T)
founderNames <- fread(founders, header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)

# ============================================================================
# Normalize the collapsed founders across replicates within treatments using cutree as the key 

tictoc::tic('total')
chem <- gsub("\\[", "", chem)
chem <- gsub("\\]", "", chem)
print(chem)
flush.console()
hapFreqsFilt <- hapFreqs[poolroot %like% chem | poolroot %like% "BAS02"]
txDT <- hapFreqsFilt[, 'id' := NULL
][, c('cp', 'id') := .(paste(chr, pos, sep = "_"), poolroot)]
chemId <- hapFreqsFilt[poolroot %like% chem, "poolroot"][1]
chem <- paste(strsplit(chemId$poolroot, "")[[1]][4:10], collapse = "")
txDT$Chemical <- treatKeyDT[V1 %like% chem, Chemical]
txDT$Replicate <- gsub(".+?(?=R[0-9][0-9])", "", txDT$id, perl = T)
txDT$Week <- paste(strsplit(unique(txDT$id)[2], "")[[1]][4:5], collapse = "")
txDT$chemWeek <- paste(txDT$Chemical, txDT$Week, sep = "_")
txDT <- txDT[, poolroot := NULL][chr != 17]
print(txDT$chemWeek[1])
flush.console()
chemLong <- dcast(txDT, chr + pos + gp ~ id, value.var = "cutree")
cutreeLists <- apply(chemLong[,4:ncol(chemLong)],1,function(x) strsplit(x, ";"))
newTrees <- lapply(cutreeLists, function(x) {
    toNums <- lapply(x, as.numeric)
    newTrees <- lapply(1:length(toNums), function(y) {
        cuttree <- toNums[[y]]
        newTree <- as.character(unlist(lapply(tapply(LETTERS[1:17],cuttree,function(x) paste(x)),function(x) paste(x,collapse=''))))
        newTree
    } )
    names(newTrees) <- names(toNums)
    newTrees
} )
print("About to merge!")
flush.console()
counter <- 0
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(newTrees), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")
mergeCutrees <- Map(merge.master, newTrees)
close(pb)
print("Merge complete")
keyCutrees <- lapply(mergeCutrees, function(x) {
    cuttree = rep(0,17)
    for(i in 1:length(x)){
        hits=unlist(strsplit(x[i],''))
        mhits=match(hits,LETTERS[1:17])
        cuttree[mhits] <- i
    }
    paste(cuttree, collapse = ";") 
} )
replaceCutrees <- do.call(rbind, keyCutrees)

txDT[, cutree := rep(replaceCutrees, uniqueN(txDT$id))]
sepCuts <- lapply(txDT$cutree, function(x) as.numeric(strsplit(x, ";", fixed = TRUE)[[1]]))
sepFreqs <- lapply(txDT$founderfreqs, function(x) as.numeric(strsplit(x, ";", fixed = TRUE)[[1]]))
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
txDT2 <- cbind(txDT, newFounderCols)
txDT2 <- txDT2[, 'founderfreqs' := NULL]
txDT2$heterozygosity <- as.numeric(txDT2$heterozygosity)
txDT2$sumFreqs <- as.numeric(txDT2$sumFreqs)
setcolorder(txDT2, c('chr', 'pos', 'gp', 'cp', 'Chemical', 'Replicate', 'id'))
fwrite(txDT2, file=paste0(projectDir, "/Tables/All_reps_tx_tables/", txDT$chemWeek[1], "_allreps_hap_freqs_DT.txt"), col.names = T, sep = "\t")
tictoc::toc()

# ============================================================================
# Trouble-shooting