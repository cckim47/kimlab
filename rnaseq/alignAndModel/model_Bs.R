## Use to get expression estimates, using BitSeq
## Uses output from transcriptome alignments in sam format
## Specify the aligner with the aligner code (Ka or Hs) in the command line call 

library(BitSeq)

## Get aligner information
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1){
	stop("Must provide alignment code.  Ex: Ka", call.=FALSE)
	}
alignCode = args[1]

## Specify transcriptome index
transcriptome <- "/opt/index/hs/genes.fa"

## Generate expression estimates in RPKM and counts and write to files
getExpression(paste0(alignCode,"/aligned.sam"), transcriptome, type="RPKM", log=FALSE, outPrefix=paste0(outDir,"/alignedExpressionRPKM"))
estimateExpression(paste0(outDir,"/alignedExpressionRPKM.prob"),paste0(outDir,"/alignedExpressionCounts"), outputType="counts", trInfoFile=paste0(outDir,"/alignedExpressionRPKM.tr"))
getMeanVariance(paste0(outDir,"/alignedExpressionCounts.counts"),paste0(outDir,"/alignedExpressionCounts.mean"))

## Collate results into a single table with ENSG, ENST, RPKM mean, counts mean, and write to file
rpkmExp <- read.table(paste0(outDir,"/alignedExpressionRPKM.mean"))
names(rpkmExp) <- c("rpkmMean","rpkmVar")
countsExp <- read.table(paste0(outDir,"/alignedExpressionCounts.mean"))
names(countsExp) <- c("countsMean","countsVar")
trList <- read.table(paste0(outDir,"/alignedExpressionRPKM.tr"))
names(trList) <- c("Empty","ENST","Length","CalcLength")
trExp <- cbind(trList,rpkmExp,countsExp)
geneTable <- read.table("../../geneConversionTable.txt", header=TRUE)
geneTrExp <- merge(trExp,geneTable, by.x="ENST", by.y="ENST")[,c(1,5,7,9)]
write.table(geneTrExp,paste0(outDir,"/transcriptExp.txt"), row.names=FALSE, quote=FALSE)
