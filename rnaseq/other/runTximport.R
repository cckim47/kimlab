## Script to convert transcript level expression to gene level expression
## Will work for Kallisto, Salmon, Sailfish abundance estimates
## Specify the RA/EM combination and the location of a gene/transcript conversion table in the command line call

library(tximport)

## Read command line argument
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2){
	stop("Supply one RA/EM name (ex. KaKa) and geneConversionTable location (ex. ../../geneConversionTable.txt)", call.=FALSE)
	}
subdir <- args[1]
geneConversion <- read.table(args[2], header=TRUE)

## Determine EM type
if (length(grep("*Ka$",subdir))==1) {
	emType = "kallisto"
	emFile = "abundance.tsv"
	} else if (length(grep("*Sl$", subdir))==1) {
	emType = "sailfish"
	emFile = "quant.sf"
	} else if (length(grep("*Sf$", subdir))==1) {
	emType = "salmon"
	emFile = "quant.sf"
	} else if (length(grep("*Sn4Su$", subdir))==1) {
	emType = "salmon"
	emFile = "quant.sf"
	} else if (length(grep("*Sl4Su$", subdir))==1) {
	emType = "sailfish"
	emFile = "quant.sf"
	} else if (length(grep("*Sn$", subdir))==1) {
	emType = "salmon"
	emFile = "aligned.counts/quant.sf"
	} else if (length(grep("*Sq$", subdir))==1) {
	emType = "salmon"
	emFile = "quant.sf"
	} else {
	stop("RA/EM combination not recognized", call.=FALSE)
	}
 
## Make transcript to gene map
t2g <- geneConversion[,c("ENST","ENSG")]

## Run tximport
txi <- tximport(paste0(subdir,"/",emFile), type = emType, tx2gene = t2g)

## Write outputs to files
write.table(txi$counts, paste0(subdir,"2g/genes.counts"), row.names=TRUE, col.names=FALSE, quote=FALSE)
write.table(txi$abundance, paste0(subdir,"2g/genes.tpm"), row.names=TRUE, col.names=FALSE, quote=FALSE)
write.table(txi$length, paste0(subdir,"2g/genes.length"), row.names=TRUE, col.names=FALSE, quote=FALSE)

