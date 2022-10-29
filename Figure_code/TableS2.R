library(readxl)
library(data.table)
setwd("/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/Experimental_Setups/")

file <- read_excel('Scoring_sexual_vs_asexual_SEE12.xlsx', 2, col_names = TRUE)
file <- as.data.frame(file)
file <- file[,-c(19:21)]
file <- file[-nrow(file),]

finding <- lapply(2:17, function(f) {
	sexuals <- grep('S', file[,f])
	return(c(f, sexuals))
} )
maxLen <- max(unlist(lapply(finding, length)))
findingCorrected <- lapply(finding, function(x) {
	length(x) <- maxLen
	return(x)
} )
findingDF <- as.data.frame(findingCorrected)
colnames(findingDF) <- findingDF[1,]
findingDF <- findingDF[-1,]

compareDFs <- lapply(findSexuals, function(x){
	comparing <- unlist(lapply(1:16, function(y) {
		key <- paste(colnames(x)[y], x[,y], sep = "_")
		key
	} ) )
	comparing
} )

asexfinding <- lapply(2:17, function(f) {
	asexuals <- grep('AD$|AD\\*$', file[,f])
	return(c(f, asexuals))
} )
asexmaxLen <- max(unlist(lapply(asexfinding, length)))
asexfindingCorrected <- lapply(asexfinding, function(x) {
	length(x) <- asexmaxLen
	return(x)
} )
asexfindingDF <- as.data.frame(asexfindingCorrected)
colnames(asexfindingDF) <- asexfindingDF[1,]
asexfindingDF <- asexfindingDF[-1,]


compareAsexDFs <- lapply(findAsexuals, function(x){
	comparing <- unlist(lapply(1:16, function(y) {
		key <- paste(colnames(x)[y], x[,y], sep = "_")
		key
	} ) )
	comparing
} )

findDH <- lapply(files, function(df) {
	finding <- lapply(2:17, function(f) {
		dhaps <- grep('H', df[,f])
		return(c(f, dhaps))
	} )
	maxLen <- max(unlist(lapply(finding, length)))
	findingCorrected <- lapply(finding, function(x) {
		length(x) <- maxLen
		return(x)
	} )
	findingDF <- as.data.frame(do.call(cbind,findingCorrected))
	colnames(findingDF) <- findingDF[1,]
	findingDF <- findingDF[-1,]
	findingDF
} )

compareDHDFs <- lapply(findDH, function(x){
	comparing <- unlist(lapply(1:16, function(y) {
		key <- paste(colnames(x)[y], x[,y], sep = "_")
		key
	} ) )
	comparing
} )

allSexual <- unlist(lapply(compareDFs, function(x) length(x[-grep('NA', x)])))
allAsexual <- unlist(lapply(compareAsexDFs, function(x) length(x[-grep('NA', x)])))
allDH <- unlist(lapply(compareDHDFs, function(x) length(x[-grep('NA', x)])))

findSexualSims <- Reduce(intersect, compareDFs)
noNasSexualSims <- findSexualSims[-grep('NA', findSexualSims)]
sexualSims <- length(noNasSexualSims)
findAsexualSims <- Reduce(intersect, compareAsexDFs)
noNasAsexualSims <- findAsexualSims[-grep('NA', findAsexualSims)]
asexualSims <- length(noNasAsexualSims)
findDHSims <- Reduce(intersect, compareDHDFs)
noNasDHSims <- findDHSims[-grep('NA', findDHSims)]
dhSims <- length(noNasDHSims)

allSamples1 <- 95+146+14
allSamples2 <- 56+170+14
allSamples3 <- 75+152+13
allSamplesTotal <- 240

fracAgreed <- (sexualSims + asexualSims + dhSims)/240
fracDH <- 12/14
fracSexual <- 

findSims12 <- intersect(compareDFs[[1]], compareDFs[[2]])
noNas12 <- findSims12[-grep('NA', findSims12)]
findSims13 <- intersect(compareDFs[[1]], compareDFs[[3]])
noNas13 <- findSims13[-grep('NA', findSims13)]
findSims23 <- intersect(compareDFs[[2]], compareDFs[[3]])
noNas23 <- findSims23[-grep('NA', findSims23)]

