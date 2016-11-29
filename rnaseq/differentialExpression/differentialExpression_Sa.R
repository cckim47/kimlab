## Rscript to run samSeq on all counts, FPKM, and TPM expression matrices
## Expression matrices have samples in columns and genes/transcripts in rows

library(samr)
library(plyr)


########################
#### Wrapper  Input ####
########################

## Set up input and output files and directories

mainDir<-getwd()
mainDirOut<-"de/differentialExpression"
sigDirOut<-"de"

conversionTable<-read.table("conversion_tgs.txt", header=TRUE)

expTypes <- c("countsGn","countsTx","fpkmGn", "fpkmTx","tpmGn", "tpmTx")

for (expType in expTypes){
  
  subDirMatrices<-paste0("em/",expType)
  file.names <- dir(file.path(mainDir,subDirMatrices), pattern =".txt")
  
  for (fileName in file.names){
    RaEm<-unlist(strsplit(fileName, split='_', fixed=TRUE))[1]
    RaEmDe<-paste0(RaEm,"Sa")
    
    ## Define output directory and create if necessary
    subDirOut<-paste0(RaEmDe,"_",expType)
    if(file.exists(file.path(mainDirOut,subDirOut))){print(paste0("Combination ", RaEmDe,expType, " has been run already. Delete or rename ALL output folders for this comparison to rerun."))
    } else{
      print(paste0("Analyzing ", RaEmDe, "_", expType, "\n"))
      dir.create(file.path(mainDirOut,subDirOut))
      
      ## Load in expression matrix 
      expMat<-read.table(file.path(mainDir,subDirMatrices,fileName), header=TRUE, row.names=1)
      


    ## Set up expression matrix and Groupings ##
    expMat[]<-as.integer(round(as.matrix(expMat)))
    groupings<-colnames(expMat)
    groupings<-as.factor(substr(groupings,1,nchar(groupings)-2))
    group1<-levels(groupings)[1]
    group2<-levels(groupings)[2]
    comparison<-paste0(group1,"v",group2)

    ## Define if comparing gene or transcript
    identifier<-substr(row.names(expMat)[1], 1, 4)
    if(identifier == "ENST"){ compType <- "transcript"}
    if(identifier == "ENSG"){ compType <- "gene"}

    ######################
    #### Run SamSeq ####
    ######################


    ## Run Sam ##
    samfit <- SAMseq(expMat, groupings, resp.type = "Two class unpaired", geneid=row.names(expMat), fdr=0.05)
       
    ## Write main output
    write.table(samfit$siggenes.table$genes.lo, sep="\t", file.path(mainDirOut,subDirOut,paste0(RaEmDe,code,"_","upGenes.txt")), quote=FALSE)
    write.table(samfit$siggenes.table$genes.up,sep="\t", file.path(mainDirOut,subDirOut,paste0(RaEmDe,code,"_","downGenes.txt")),quote=FALSE)
       
    ## Reformat Output for precision and Recall 
       
    samUp<-read.table(file.path(mainDirOut,subDirOut,paste0(RaEmDe,code,"_","upGenes.txt")), header=TRUE, sep="\t")
    samDown<-read.table(file.path(mainDirOut,subDirOut,paste0(RaEmDe,code,"_","downGenes.txt")), header=TRUE, sep="\t")
    allSam<-rbind(samUp,samDown)
    names(allSam)[names(allSam)=="Gene.Name"] <- identifier
    allSam<-rename(allSam, c("Gene.Name" = "UNIQID"))
    allSam<-merge(allSam,conversionTable,by= identifier)
       

    if(compType == "gene"){
	write.table(unique(allSam), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigGenes05.txt")),quote=FALSE)
	} else {
	    write.table(unique(allSam), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigTranscripts05.txt")),quote=FALSE)
	  }
    write.table(unique(allSam$geneName),sep="\t", row.names=FALSE,file.path(sigDirOut,expType,paste0(RaEmDe,"_","sigSymbol.txt")), quote=FALSE, col.names=FALSE)

    }
  }
}

