#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Annotate mutations
# Description: Further annotate the list of potential mutations using KEGG/GO terms

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: analysis_step15c_annotate_mutations_KEGG_GO.R <mutations_ann> <helper>", call.=FALSE)
}

mutations_ann <-  args[1]
helper <- args[2]

# ============================================================================
# Load packages and sourced files

library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(RCurl)
library(tictoc)
#library(RIPSeeker)
library(dplyr)
library(data.table)
library(org.Sc.sgd.db)
library(GO.db)
#library(KEGG.db)
library(KEGGREST)
library(pathview)
source(helper)
#source("/Users/robertlinder/Dropbox/Long_lab/SEE01/Linder-et-al-2022/bin/helper_functions.R")
# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

# ============================================================================
# Load data

mutsDT <- fread(mutations_ann)
#mutsDT <- fread("/Users/robertlinder/Dropbox/Long_lab/SEE01/Linder-et-al-2022/results/All_mutations_annotated_less_stringent2.txt")

# ============================================================================
# Adding GO and KEGG PATH terms to genes

mutsDTypd <- mutsDT[chemical == "YP"]
mutsDTtreat <- mutsDT[!mutID %in% unique(mutsDTypd$mutID)]

organism <- "saccharomyces cerevisiae"
matches <- unlist(sapply(1:ncol(korg), function(i) {agrep(organism, korg[, i])}))
mutsDT$GO <- NA
mutsDT$PATH <- NA

counter <- 0
addGOterms <- lapply(mutsDT$gene_hit, function(gene) {
    counter <<- counter + 1
    print(counter)
    flush.console()
    multSplit <- unlist(strsplit(gene, ";")[[1]])
    geneList <- unlist(lapply(multSplit, function(x) {
        if(grepl("-", x)) {
            location <- unlist(gregexpr("-", x))
            if(length(location) > 1) {
                firstSplit <- strsplit(x, "-")[[1]]
                firstPaste <- paste0(firstSplit[1], "-", firstSplit[2])
                splitGenes <- c(firstPaste, firstSplit[3])
            }
            else if((nchar(x) - location) > 2) {
                splitGenes <- strsplit(x, "-")[[1]]
            } else {
                splitGenes <- x
            }
        } else {
            splitGenes <- x
        }
    return(splitGenes)
    } ) )
    #genesHit <- which(unlist(lapply(charList, nchar)) > 2)
    #geneList <- charList[genesHit]
    orfKeys <- geneList[geneList %in% keys(org.Sc.sgd.db)]
    commonKeys <- geneList[geneList %in% keys(org.Sc.sgd.db, keytype = "COMMON")]
    if(length(orfKeys) > 0 & length(commonKeys) > 0) {
        annotOrfs <- select(org.Sc.sgd.db, keys = orfKeys, columns = c("GO", "ORF"))
        annotCommon <- select(org.Sc.sgd.db, keys = commonKeys, columns = c("GO", "ORF"), keytype = "COMMON")
        names(annotCommon)[1] <- "GENE"
        names(annotOrfs)[1] <- "GENE"
        if(!"ORF" %in% names(annotOrfs)) {annotOrfs$ORF <- NA}
        allAnnotGenes <- rbind(annotOrfs, annotCommon) } else if(length(orfKeys) > 0) {
            allAnnotGenes <- select(org.Sc.sgd.db, keys = orfKeys, columns = c("GO", "ORF"))
            names(allAnnotGenes)[1] <- "GENE"
            if(!"ORF" %in% names(allAnnotGenes)) {allAnnotGenes$ORF <- NA}} else {
                allAnnotGenes <- select(org.Sc.sgd.db, keys = commonKeys, columns = c("GO", "ORF"), keytype = "COMMON")
                names(allAnnotGenes)[1] <- "GENE"
                if(!"ORF" %in% names(allAnnotGenes)) {allAnnotGenes$ORF <- NA}}
    findNAs <- which(is.na(allAnnotGenes$SGD))
    removeNAs <- unlist(lapply(findNAs, function(rm) {
        if(length(which(allAnnotGenes$ORF == allAnnotGenes$ORF[rm])) > 1) {return(rm)} 
    } ) )
    if(length(removeNAs) > 0) {allAnnotGenes <- allAnnotGenes[-removeNAs,]} else{
        allAnnotGenes <- allAnnotGenes}
    #allAnnotGenes$PATH <- unlist(lapply(allAnnotGenes$PATH, function(x) {paste0('sce', x)}))
    findGOs <- allAnnotGenes$GO[!(is.na(allAnnotGenes$GO))] 
    if(length(findGOs > 0)) {
        annotGOs <- select(GO.db, keys = allAnnotGenes$GO, columns = "TERM",  keytype = "GOID")
        findDups <- which(duplicated(annotGOs$GOID))
        if(length(findDups) > 0 ){goIds <- annotGOs[-findDups,]} else{goIds <- annotGOs}
        goIds <- goIds[!(is.na(goIds$GOID)),]
        allAnnotGenes$GO[allAnnotGenes$EVIDENCE == "ND"] <- NA
        allAnnotGenes$GO[!is.na(allAnnotGenes$GO)] <- unlist(lapply(allAnnotGenes$GO[!is.na(allAnnotGenes$GO)], function(go) {
            goTerm <- goIds$TERM[goIds$GOID == go]
            if(length(goTerm) > 1) {return(goTerm[1])} else{
                return(goTerm)}
        } ) )}
    if(nrow(allAnnotGenes) == 0) {
        allAnnotGenes[1,] <- NA
        allAnnotGenes$GENE <- geneList}
    dubious <- which(is.na(allAnnotGenes$ORF))
    nonDubious <- which(!is.na(allAnnotGenes$ORF))
    if(length(dubious) > 0) {allNullPaths <- data.frame(gene = allAnnotGenes$GENE[dubious], paths = NA, stringsAsFactors = FALSE)}
    findKeggId <- names(unlist(lapply(unique(allAnnotGenes$ORF[!(is.na(allAnnotGenes$ORF))]), keggFind, database = "sce")))
    #annotPaths <- unique(allAnnotGenes$PATH[!grepl("^sceNA$", allAnnotGenes$PATH)])
    if(length(findKeggId) > 0) { ### need to modify so deals with cases where are NAs for dubious ORFs and values for other genes
        if(length(findKeggId) <= 10) {
            addPath <- keggGet(findKeggId)
            extractPath <- lapply(addPath, function(path) {
                name <- which(names(path) == "ENTRY")
                name2 <- which(names(path) == "NAME")
                if(length(name2) < 1) {name2 <- NA}
                if(any(allAnnotGenes$GENE %in% path[[name2]]) | any(allAnnotGenes$ORF %in% path[[name]])) {
                    findPath <- which(names(path) == "PATHWAY")
                    findBrite <- which(names(path) == "BRITE")
                    findDefinition <- which(names(path) == "DEFINITION")
                    if(length(findPath) > 0) {
                        map <- path[findPath]} else if (length(findBrite) > 0 & length(findPath) == 0) {
                            map <- path[[findBrite]][4]
                            map <- strsplit(map, "[0-9]|\\[")[[1]][6] } else if(length(findDefinition) > 0 & length(findPath) == 0 & length(findBrite) == 0) {
                                map <- path[[findDefinition]]} else {
                                    map <- NA}
                    entries <- unlist(c(path[name], map))
                    data.frame(gene = entries[1], paths = entries[2:length(entries)], stringsAsFactors = FALSE) } } ) 
            extractPath <- Filter(Negate(function(i) is.null(unlist(i))), extractPath)
            allPaths <- do.call(rbind, extractPath)} else {
                splitPaths <- split(findKeggId, (seq(length(findKeggId))-1) %/% 10)
                splitLoops <- lapply(splitPaths, function(sp) {
                    addPath <- keggGet(sp)
                    pathCounter <- 0
                    extractPath <- lapply(addPath, function(path) {
                        pathCounter <<- pathCounter + 1
                        print(pathCounter)
                        flush.console()
                        name <- which(names(path) == "ENTRY")
                        name2 <- which(names(path) == "NAME")
                        if(length(name2) < 1) {name2 <- NA}
                        if(any(allAnnotGenes$GENE %in% path[[name2]]) | any(allAnnotGenes$ORF %in% path[[name]])) {
                            findPath <- which(names(path) == "PATHWAY")
                            findBrite <- which(names(path) == "BRITE")
                            findDefinition <- which(names(path) == "DEFINITION")
                            if(length(findPath) > 0) {
                                map <- path[findPath]} else if(length(findBrite) > 0 & length(findPath) == 0) {
                                    map <- path[[findBrite]][4]
                                    map <- strsplit(map, "[0-9]|\\[")[[1]][6]} else if(length(findDefinition) > 0 & length(findPath) == 0 & length(findBrite) == 0) {
                                        map <- path[[findDefinition]]} else{
                                            map <- NA}
                            entries <- unlist(c(path[name], map))
                            data.frame(gene = entries[1], paths = entries[2:length(entries)], stringsAsFactors = FALSE) } } )
                    extractPath <- Filter(Negate(function(i) is.null(unlist(i))), extractPath)
                    do.call(rbind, extractPath) } ) 
                allPaths <- do.call(rbind, splitLoops) } } else if(!exists("allNullPaths")){
                    allNullPaths <- data.frame(gene = unique(allAnnotGenes$GENE), paths = NA, stringsAsFactors = FALSE)}
    if(exists("allNullPaths") & exists("allPaths")) {
        allPaths <- rbind(allNullPaths, allPaths)} else if(exists("allNullPaths")) {
            allPaths <- allNullPaths} else if(exists("allPaths")) {
                allPaths <- allPaths}
    #allAnnotGenes$PATH[grepl("^sceNA$", allAnnotGenes$PATH)] <- NA
    addingBackGenes <- geneList[!geneList %in% allAnnotGenes$GENE] ## fix this part
    if(length(addingBackGenes) > 0 & !exists("allPaths")){
        allAnnotGenes[(nrow(allAnnotGenes)+1):(nrow(allAnnotGenes)+length(addingBackGenes)),1] <- addingBackGenes
        allPaths <- data.frame(gene = allAnnotGenes$GENE, paths = NA, stringsAsFactors = FALSE)} else if(length(addingBackGenes) > 0 & exists("allPaths")) {
            allAnnotGenes[(nrow(allAnnotGenes)+1):(nrow(allAnnotGenes)+length(addingBackGenes)),1] <- addingBackGenes} ### have to add all missing ones back here
    combineGOTerms <- split(allAnnotGenes, allAnnotGenes$GENE)
    reorder <- match(geneList, names(combineGOTerms))
    allAnnotGenes$ORF[is.na(allAnnotGenes$ORF)] <- allAnnotGenes$GENE[is.na(allAnnotGenes$ORF)]
    combinePATHTerms <- split(allPaths, allPaths$gene)
    reorderPath <- match(unique(allAnnotGenes$ORF), names(combinePATHTerms))
    combineGOTermsOrdered <- combineGOTerms[c(reorder)]
    combinePathTermsOrdered <- combinePATHTerms[c(reorderPath)]
    findUniqueGO <- paste(unlist(lapply(combineGOTermsOrdered, function(g) {
        anyNAs <- which(!is.na(unique(g$GO)))
        if(length(anyNAs) > 0) {return(unique(g$GO)[anyNAs[1]])} else {
            return(unique(g$GO[1]))}})), collapse = "; ") ## only listing one GO term per gene
    findUniquePATH <- paste(unlist(lapply(combinePathTermsOrdered, function(g) {
        anyNAs <- which(!is.na(unique(g$paths)))
        if(length(anyNAs) > 0) {return(unique(g$paths)[anyNAs[1]])} else {
            return(unique(g$paths[1]))}})), collapse = "; ")
    mutsDT$GO[counter] <- findUniqueGO
    mutsDT$PATH[counter] <- findUniquePATH
    mutsDT[counter,]
} )

mutsDTcomp <- do.call(rbind, addGOterms)
mutsDTcomp <- as.data.table(mutsDTcomp)
fwrite(mutsDTcomp, "SNVsAnnotatedAbr.txt", quote = F, sep = "\t", col.names = T)
# ============================================================================
# Trouble-shooting