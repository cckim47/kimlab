## Script to convert transcript level expression to gene level expression
## Will work for BitSeq abundance estimates only
## Must provide the RA/EM directory name, as well as the location of a transcript/gene conversion table in the command line script call
## Uses same method as tximport, but variables must be set up manually, as bitseq format is not accepted by tximport

library(tximport)

## Read command line argument
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2){
	stop("Supply one RA/EM name (ex. KaBs) and geneConversionTable location (ex. ../../geneConversionTable.txt)", call.=FALSE)
	}
subdir <- args[1]
geneConversion <- read.table(args[2], header=TRUE)

## Check EM type
if (length(grep("*Bs$",subdir))!=1) {
	stop("RA/EM combination not recognized", call.=FALSE)
	}
 
## Make transcript to gene map
t2g <- geneConversion[,c("ENST","ENSG")]

## Import data, format as done in tximport, and run tximport
transcriptExp = read.table(paste0(subdir,"/transcriptExp.txt"), header = TRUE)
lengthTable = read.table(paste0(subdir,"/alignedExpressionRPKM.tr"), header = FALSE)
names(lengthTable) = c("Blank", "ENST", "Length", "EffLength")
transExpLen = merge(transcriptExp, lengthTable, by = "ENST")
mat <- matrix(nrow = nrow(transExpLen), ncol = 1)
rownames(mat) = transExpLen[["ENST"]]
colnames(mat) = "transExpName"
abundanceMatTx <- mat
countsMatTx <- mat
lengthMatTx <- mat
abundanceMatTx[,1] <- transExpLen[["rpkmMean"]]
countsMatTx[,1] <- transExpLen[["countsMean"]]
lengthMatTx[,1] <- transExpLen[["Length"]]
txi <- list(abundance = abundanceMatTx, counts = countsMatTx, length = lengthMatTx, countsFromAbundance = "no")
txiGene <- summarizeToGene(txi, tx2gen = t2g)


## Write outputs to files
write.table(txiGene$counts, paste0(subdir,"2g/genes.counts"), row.names=TRUE, col.names=FALSE, quote=FALSE)
write.table(txiGene$abundance, paste0(subdir,"2g/genes.rpkm"), row.names=TRUE, col.names=FALSE, quote=FALSE)
write.table(txiGene$length, paste0(subdir,"2g/genes.length"), row.names=TRUE, col.names=FALSE, quote=FALSE)

