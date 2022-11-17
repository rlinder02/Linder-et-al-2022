#!/usr/bin/env Rscript

##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2022
# Title:       Annotate mutations
# Description: Further annotate the list of potential mutations using the S. cerevisiae gff file and motifbreakr

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    stop("Usage: analysis_step15b_annotate_mutations_post_SNPeff.R <ann_mutations> <gff> <helper> <old_muts> <new_muts> <offsets>", call.=FALSE)
}

ann_mutations <-  args[1]
gff <- args[2]
helper <- args[3]
old_muts <-  args[4]
new_muts <- args[5]
offsets <- args[6]

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
install.packages("tictoc", repos = "http://cran.us.r-project.org")

# ============================================================================
# Load packages and sourced files

library(GenomicRanges)
library(rtracklayer)
library(motifbreakR)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(MotifDb)
library(RCurl)
library(tictoc)
#library(RIPSeeker)
#library(dplyr)
library(data.table)
#library(org.Sc.sgd.db)
#library(pathview)
source(helper)

# ============================================================================
# Set global options

projectDir <- getwd()

# ============================================================================
# Custom functions

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
# Load data

offsets <- fread(offsets, header = T)
g_l <- posoff_chr.bounds(offsets)
mutsAnnDF <- read.delim(ann_mutations, skip = 12, header = F, check.names = FALSE, colClasses = c("character", "integer", "character", "character", "character", "character", "character", "character"))
oldMuts <- fread(old_muts, header = T)
rhcMuts <- oldMuts[sample %like% "^RHC"]
baseMuts <- oldMuts[sample %like% "^postsp"]
newMuts <- fread(new_muts, header = T)
newMuts <- rbind(newMuts, rhcMuts, baseMuts)


# ============================================================================
# Converting gff files to granges object

tags_to_keep <- c("X_element", "telomeric_repeat", "gene", "ARS", "long_terminal_repeat", "ncRNA_gene", "tRNA_gene", "snoRNA_gene", "centromere", "LTR_retrotransposon", "transposable_element_gene", "pseudogene", "Y_prime_element", "telomerase_RNA_gene", "snRNA_gene", "silent_mating_type_cassette_array", "mating_type_region", "intein_encoding_region", "rRNA_gene", "external_transcribed_spacer_region", "internal_transcribed_spacer_region", "non_transcribed_region", "origin_of_replication")
myGranges <- gff_to_granges(gff, tags_to_keep)

# ============================================================================
# Need to rerun this to regenerate the treatMutsCp table 

newMuts <- newMuts[,chr2 := as.integer(as.roman(substr(chr, 4, nchar(chr))))
][chr2 == 1000, chr2 := 17
][order(chr2, pos, sample)]
newMuts <- newMuts[, gp := .(pos + g_l[chr2])]
newMuts <- newMuts[, id := sample][, sample := NULL]
newMuts$chemical <- unlist(lapply(1:nrow(newMuts), function(x){
    chemName <- paste(strsplit(newMuts$id[x], "")[[1]][9:10], collapse = "")
    return(chemName)
} ) )
newMuts <- newMuts[,mutID := paste(chr, pos, alt, sep = "_")]
newMuts <- newMuts[, c("fraction", "totalCov") := .(Nalt/(Nalt+Nref), Nref+Nalt)]
firstFilter <- newMuts[id %like% "^A[0-9]" | id %like% "^B[0-9]" | id %like% "^AB[0-9]" | id %like% "^BAS02$" | id %like% "^RHC" | id %like% "^DIP02$" | id %like% "^Base" | id %like% "FA02P" | id %like% "FA11P" | id %like% "^postsp"] # Filter out all mutations that show up in the original founders (including on plate controls - FA02 and FA11 founders here) as well as in any of the base populations, including clones directly derived from the base populations, as these mutations are likely misbehaving SNPs or were missed when sequencing the original founders or occur as part of the normal culturing of cells and sample prep for sequencing; are most interested in mutations that occur for specific treatments
newMuts <- newMuts[!mutID %in% unique(firstFilter$mutID)]
newMuts <- newMuts[totalCov >= 10] # Need at least 10x coverage at the mutations' position
newMuts <- newMuts[fraction >= 0.2] # At least 20% of the calls have to be for the mutation
countMuts <- newMuts[, .N, by = c("mutID", "chemical")]
uniqueMuts <- countMuts[, uniqueN(chemical), by = mutID]
goodMuts <- uniqueMuts[V1 <= 5] # Mutations that show up in more than five different treatments are unlikely to be specific for the treatment- more likely an artifact of the culturing conditions
goodMutsDT <- newMuts[mutID %in% goodMuts$mutID]
treatMutsCp <- copy(goodMutsDT)
treatMutsCp[chr == 'chrM', chr := 'chrmt']
treatMutsCp[, AllID := paste(chr, pos, id, ref, alt, sep = "_")]

# ============================================================================
# Filter and reformat snpeff annotated mutations table 

mutsAnnDT <- as.data.table(mutsAnnDF)
setnames(mutsAnnDT, c("Chr", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
mutsAnnDT[Chr == 'chrMito', Chr := 'chrmt']
mutsAnnDT[, AllID := paste(Chr, POS, ID, REF, ALT, sep = "_")]
codingMutIdx <- "synonymous_variant|missense_variant|intron_variant|stop_gained|start_lost|stop_lost|stop_lost&splice_region_variant|splice_region_variant&stop_retained_variant|missense_variant&splice_region_variant|splice_region_variant&intron_variant|initiator_codon_variant"
noncodingMutIdx <- "downstream_gene_variant|intergenic_region|upstream_gene_variant|non_coding_transcript_exon_variant"
intergenicDT <- mutsAnnDT[INFO %like% "intergenic_region"]
downstreamDT <- mutsAnnDT[INFO %like% "downstream_gene_variant"]
upstreamDT <- mutsAnnDT[INFO %like% "upstream_gene_variant"]
codingDT <- mutsAnnDT[INFO %like% codingMutIdx]
noncodingDT <- mutsAnnDT[INFO %like% noncodingMutIdx]
mutsAnnDT2 <- mutsAnnDT[!(AllID %in% codingDT$AllID & AllID %in% noncodingDT$AllID & INFO %like% noncodingMutIdx)]
mutsAnnDT3 <- mutsAnnDT2[!(AllID %in% intergenicDT$AllID & AllID %in% downstreamDT$AllID & INFO %like% "intergenic_region")]
mutsAnnDT4 <- mutsAnnDT3[!(AllID %in% intergenicDT$AllID & AllID %in% upstreamDT$AllID & INFO %like% "intergenic_region")]
mutsAnnDT4$INFO <- unlist(lapply(mutsAnnDT4$INFO, function(x) {
    if(grepl("^ANN", x)) {return(x)} else {
        newAnn <- strsplit(x, ";")[[1]]
        newAnn2 <- paste(newAnn[2:length(newAnn)], sep = ";")
        return(newAnn2)}
} ))
codingMutsDT <- mutsAnnDT4[INFO %like% codingMutIdx]
noncodingMutsDT <- mutsAnnDT4[INFO %like% noncodingMutIdx]

# ============================================================================
# Adding annotations from S288C gff to noncoding mutations table 

noncodingMutDF <- as.data.frame(noncodingMutsDT)
noncodingMutDF$chrpos <- paste(noncodingMutDF$Chr,noncodingMutDF$POS, sep = "_")
chrposUnique <- unique(noncodingMutDF$chrpos)
chrUnique <- unlist(lapply(chrposUnique, function(x) {strsplit(x, "_")[[1]][1]}))
posUnique <- unlist(lapply(chrposUnique, function(x) {strsplit(x, "_")[[1]][2]}))
coordsDF <- data.frame(chrom = chrUnique, start = posUnique, end = posUnique, stringsAsFactors = FALSE)
coordsDF$chrpos <- paste(coordsDF$chrom, coordsDF$start, sep = "_")
coordsGR <- makeGRangesFromDataFrame(coordsDF, keep.extra.columns = TRUE)
overlappingRegion <- subsetByOverlaps(myGranges, coordsGR)
overlappingRegion2 <- subsetByOverlaps(coordsGR,myGranges)
hits <- findOverlaps(myGranges, coordsGR)
hits <- data.frame(hits)
all_hits <- lapply(unique(hits$queryHits), function(x) {
    find_hits <- hits$subjectHits[hits$queryHits == x]
    print(find_hits)
    flush.console()
    keep_values <- coordsDF$chrpos[find_hits]
    keep_values
})
allHitsDF <- data.frame(Idx = I(all_hits))
mcols(overlappingRegion) <- c(mcols(overlappingRegion), allHitsDF)
noncodingHitsDF <- data.frame(symbol = overlappingRegion$symbol, Idx = overlappingRegion$Idx, stringsAsFactors = FALSE)
findIdx <- lapply(unique(noncodingMutDF$chrpos), function(x) {
    noncodinghitsidx <- grep(x, noncodingHitsDF[,2])
    if(length(noncodinghitsidx > 0)) {
        noncodingmutationidx <- noncodingMutDF[noncodingMutDF$chrpos == x,]
        find_dups <- which(duplicated(noncodingmutationidx$AllID))
        if(length(find_dups) > 0) {noncodingmutationidx <- noncodingmutationidx[-find_dups,] } else {
            noncodingmutationidx <- noncodingmutationidx
        }
        noncodingmutationidx$INFO <- lapply(noncodingmutationidx$INFO, function(y) {
            split_info <- strsplit(y, "\\|")[[1]]
            if(length(noncodinghitsidx) > 1) {split_info[4] <- paste(noncodingHitsDF[,1][c(noncodinghitsidx)], collapse = ";")} else {
                split_info[4] <- noncodingHitsDF[,1][noncodinghitsidx]}
            split_info[c(2,5:15)] <- ""
            split_info <- paste(split_info, collapse = "|")
            split_info
        } )
        noncoding_idx_df <- do.call(rbind, noncodingmutationidx$INFO)
        noncoding_idx_df <- data.frame(INFO = noncoding_idx_df, stringsAsFactors = FALSE)
        allid <- data.frame(AllID = noncodingmutationidx$AllID, stringsAsFactors = FALSE)
        noncoding_idx_all_df <- cbind(noncoding_idx_df, allid)
        noncoding_idx_all_df
    }
} )
find_idx <- Filter(Negate(is.null), findIdx)
find_idx_df <- do.call(rbind, find_idx)
mutsAnnDT4$INFO <- unlist(lapply(1:nrow(mutsAnnDT4), function(x) {
    if(mutsAnnDT4$AllID[x] %in% find_idx_df$AllID) {
        print(x)
        new_info_idx <- which(find_idx_df$AllID == mutsAnnDT4$AllID[x])
        return(find_idx_df$INFO[new_info_idx]) } else {
            return(mutsAnnDT4$INFO[x])}
} ) )
mutsAnnDT4$Chr[mutsAnnDT4$Chr == 'chrmt'] <- 'chrM'
print("about to use motifbreakr")
flush.console()
# ============================================================================
# Adding annotations using motifbreakr for TFBS disruption for the remaining noncoding mutations

notAnnDT <- mutsAnnDT4[INFO %like% noncodingMutIdx]
notAnnDT$start <- notAnnDT$POS - 1
notAnnDT$end <- notAnnDT$POS 
notAnnDT$name <- paste(notAnnDT$Chr, notAnnDT$POS, notAnnDT$REF, notAnnDT$ALT, sep = ":")
notAnnDT$score <- 0
notAnnDT$strand <- '+'
tfbsDF <- data.frame(chromosome = notAnnDT$Chr, start = notAnnDT$start, end = notAnnDT$end, name = notAnnDT$name, score = notAnnDT$score, strand = notAnnDT$strand, stringsAsFactors = FALSE)
names(tfbsDF) <- NULL
write.table(tfbsDF, "/usr/local/lib/R/library/motifbreakR/extdata/notann2.bed", quote = F, sep = '\t', col.names = F, row.names = F)
muts_bedfile_cer <- system.file("extdata", "notann2.bed", package = "motifbreakR")
muts_from_bed_df <- read.table(muts_bedfile_cer, header = FALSE)
factors <- sapply(muts_from_bed_df, is.factor)
muts_from_bed_df[factors] <- lapply(muts_from_bed_df[factors], as.character)
muts_from_bed <- snps.from.file(file = muts_bedfile_cer, search.genome = BSgenome.Scerevisiae.UCSC.sacCer3, format = "bed")
motifs <- query(MotifDb, 'cerevisiae')
print("Got to here!")
flush.console()
tic()
motifs_disrupted <- motifbreakR(snpList = muts_from_bed, filterp = TRUE, pwmList = motifs, threshold = 1e-4, method = "ic", bkg = c(A=0.25, C=0.25, G=0.25, T=0.25), BPPARAM = BiocParallel::MulticoreParam())
toc()
motifs_disrupted_df <- data.frame(motifs_disrupted, stringsAsFactors = FALSE)
motifs_disrupted_df$seqnames <- as.character(motifs_disrupted_df$seqnames)
motifs_disrupted_df$strand <- as.character(motifs_disrupted_df$strand)
rownames(motifs_disrupted_df) <- as.numeric(1:nrow(motifs_disrupted_df))
not_ann_df <- as.data.frame(notAnnDT)
#dup_rows <- which(duplicated(rownames(motifs_disrupted_df)))
#write.table(motifs_disrupted_df, "TFBS_disrupted_noncoding_muts_df.txt", quote = F, sep = "\t")
#motifs_disrupted_df <- read.delim("TFBS_disrupted_noncoding_muts_df.txt")
#i <- sapply(motifs_disrupted_df, is.factor)
#motifs_disrupted_df[i] <- lapply(motifs_disrupted_df[i], as.character)
motifs_disrupted_df <- motifs_disrupted_df[motifs_disrupted_df$effect == "strong",]
motifs_disrupted_df$name <- motifs_disrupted_df$SNP_id
find_indices <- lapply(unique(notAnnDT$name), function(x) {
    noncodinghitsidx <- grep(x, motifs_disrupted_df$name)
    if(length(noncodinghitsidx > 0)) {
        noncodingmutationidx <- not_ann_df[not_ann_df$name == x,]
        noncodingmutationidx$INFO <- lapply(noncodingmutationidx$INFO, function(y) {
            split_info <- strsplit(y, "\\|")[[1]]
            if(length(noncodinghitsidx) > 1) {split_info[5] <- paste(unique(motifs_disrupted_df$geneSymbol[c(noncodinghitsidx)]), collapse = ";")} else {
                split_info[5] <- motifs_disrupted_df$geneSymbol[noncodinghitsidx]}
            split_info[3] <- "TFBS_Disruption"
            split_info[c(6:15)] <- ""
            split_info <- paste(split_info, collapse = "|")
            #print(split_info)
            split_info
        } )
        noncoding_idx_df <- do.call(rbind, noncodingmutationidx$INFO)
        noncoding_idx_df <- data.frame(INFO = noncoding_idx_df, stringsAsFactors = FALSE)
        allid <- data.frame(AllID = noncodingmutationidx$AllID, stringsAsFactors = FALSE)
        noncoding_idx_all_df <- cbind(noncoding_idx_df, allid)
        noncoding_idx_all_df
    }
} )
find_indices <- Filter(Negate(is.null), find_indices)
find_indices_df <- do.call(rbind, find_indices)
identifier <- lapply(strsplit(find_indices_df$INFO, "\\|"), function(x) {paste(x[c(1:2,4)], collapse = "|")})
identifier <- do.call(rbind, identifier)
find_indices_df$AlllID <- paste(find_indices_df$AllID, identifier, sep = "_")
all_muts_df1 <- as.data.frame(mutsAnnDT4)
#all_muts_df1$INFO[!grepl("^ANN", all_muts_df1$INFO)] <- "ANN=A|start_lost|HIGH|MMR1|YLR190W|transcript|YLR190W|protein_coding|1/1|c.3G>A|p.Met1?|3/1476|3/1476|1/491||"
identifier2 <- lapply(strsplit(all_muts_df1$INFO, "\\|"), function(x) {paste(x[c(1:2,4)], collapse = "|")})
identifier2 <- do.call(rbind, identifier2)
all_muts_df1$AlllID <- paste(all_muts_df1$AllID, identifier2, sep = "_")
all_muts_df1$INFO <- unlist(lapply(1:nrow(all_muts_df1), function(x) {
    if(all_muts_df1$AlllID[x] %in% find_indices_df$AlllID) {
        #print(x)
        new_info_idx <- which(find_indices_df$AlllID == all_muts_df1$AlllID[x])
        #if(length(new_info_idx) > 1) {print(x)}
        return(find_indices_df$INFO[new_info_idx]) } else {
            return(all_muts_df1$INFO[x])}
} ) )
all_muts_df <- all_muts_df1[,c(1:9)]
all_muts_df$gene_hit <- unlist(lapply(all_muts_df$INFO, function(x)  {
    gene <- strsplit(x, "\\|")[[1]][4]
    gene 
} ) )
all_muts_df$tx<- unlist(lapply(all_muts_df$ID, function(x)  {
    tx <- paste(strsplit(x, "")[[1]][1:13], collapse = "")
    tx 
} ) )
all_muts_df$dose<- unlist(lapply(all_muts_df$ID, function(x)  {
    gene <- paste(strsplit(x, "")[[1]][9:13], collapse = "")
    gene 
} ) )

all_muts_df$tx_gene <- paste(all_muts_df$gene_hit, all_muts_df$tx, sep = "_")
all_muts_df$dose_gene <- paste(all_muts_df$gene_hit, all_muts_df$dose, sep = "_")
mutsDT <- as.data.table(all_muts_df)
mutsDT2 <- mutsDT[treatMutsCp[,.SD, .SDcols = c("AllID", "chemical", "fraction", "totalCov", "gp")], on = .(AllID), nomatch = 0]
mutsDT2 <- mutsDT2[order(ID, Chr, POS)]
mutsDT2$mut_type <- unlist(lapply(mutsDT2$INFO, function(x) {
    descrip <- strsplit(x, "\\|")[[1]]
    if(any(grepl("TFBS_Disruption", descrip))) {
        return("TFBS_Disruption")} else if(any(grepl("MODIFIER", descrip))) {
            return("Intergenic")} else {
                return(descrip[2])}
} ) )
mutsDT2[mut_type == "missense_variant", mut_type := "Missense"
        ][mut_type == "synonymous_variant", mut_type := "Synonymous"
          ][mut_type == "stop_lost&splice_region_variant", mut_type := "Stop_lost"
            ][mut_type == "stop_gained", mut_type := "Premature_stop"
              ][mut_type == "splice_region_variant&intron_variant", mut_type := "Splice_variant"
                ][mut_type == "start_lost", mut_type := "Start_lost"
                  ][mut_type == "splice_region_variant&stop_retained_variant", mut_type := "Splice&stop_variant"
                    ][mut_type == "stop_retained_variant", mut_type := "Stop_variant"]
mutsDT2[, mutID := paste(Chr, POS, REF, ALT, sep = "_")]
fwrite(mutsDT2, "All_mutations_annotated_less_stringent2.txt", sep = "\t", col.names = T, row.names = F, quote = F)
#all_muts_df_reduced <- mutsDT2[,c(1,2,3,4,5,10,12), with = FALSE]
#write.table(all_muts_df_reduced, "All_mutations_reduced_annotated_less_stringent2.txt", sep = "\t", col.names = T, row.names = F, quote = F)
# ============================================================================
# Trouble-shooting