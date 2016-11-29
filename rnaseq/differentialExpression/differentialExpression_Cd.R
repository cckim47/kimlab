## Rscript to run Cuffdiff on all tophat outputs

#######################
#### Wrapper Input ####
#######################

## Specify index and run information
genes <- "/opt/index/hs/genes.gtf"
genome <- "/opt/index/hs/genome.fa"
cpu <- 8

## Specify sample specific data
mainDir <- getwd()
controlGroup <- "classical"
testGroup <- "nonclassical"
alignCode <- "Th"
mainDirIn <- "mapAndModel"

mainDirOut <- "de/differentialExpression"
sigDirOut <- "de"

RaEmDe <- "ThCuCd"
expType <- "fpkmGn"

## Define output directory and create if necessary
subDirOut<-paste0(RaEmDe,"_",expType)
if(file.exists(file.path(mainDirOut,subDirOut))){print(paste0("Combination ", RaEmDe,expType, " has been run already. Delete or rename ALL output folders for this comparison to rerun."))
} else{
  print(paste0("Analyzing ", RaEmDe, "_", expType, "\n"))
  dir.create(file.path(mainDirOut,subDirOut))
  
  ## Specify location of tophat results and samples
  control_id <- dir(file.path(mainDir, mainDirIn), pattern=paste0("^",controlGroup,"*"))
  control_dirs <- sapply(control_id, function(id) file.path(mainDir, mainDirIn, id, alignCode, "accepted_hits.bam")) 
  controlString <- paste(control_dirs, collapse=",")
  test_id <- dir(file.path(mainDir, mainDirIn), pattern=paste0("^", testGroup,"*"))
  test_dirs <- sapply(test_id, function(id) file.path(mainDir, mainDirIn, id, alignCode, "accepted_hits.bam")) 
  testString <- paste(test_dirs, collapse=",")  
  
  ######################
  #### Run Cuffdiff ####
  ######################
  
  system(paste0("cuffdiff -o ", mainDirOut, "/", subDirOut, " -b ", genome, " -p ", cpu, 
                " -L ", controlGroup, ",", testGroup, " -u ", genes, " ", controlString, " ", testString))
  
  #######################
  #### Write Outputs ####
  #######################
  
  ## First, read in cuffdiff results
  genesRes <- read.table(file.path(mainDir, mainDirOut, subDirOut, "gene_exp.diff"), header=TRUE)
  transcriptRes <- read.table(file.path(mainDir, mainDirOut, subDirOut, "isoform_exp.diff"), header=TRUE)
  
  ## Identify significant genes / transcripts
  genesSig <- subset(genesRes, significant=="yes")
  transcriptsSig <- subset(transcriptRes, significant=="yes")
  
  ## Write results
  write.table(genesSig$gene, sep="\t", row.names=FALSE, file.path(sigDirOut, expType, paste0(RaEmDe, "_sigSymbol.txt")), quote=FALSE, col.names=FALSE)
  write.table(unique(transcriptsSig$gene),sep="\t", row.names=FALSE, file.path(sigDirOut, "fpkmTx", paste0(RaEmDe, "_sigSymbol.txt")), quote=FALSE, col.names=FALSE)
}
