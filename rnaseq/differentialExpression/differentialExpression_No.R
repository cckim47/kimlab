## Rscript to run NOISeqBIO on all counts, FPKM, and TPM expression matrices
## Expression matrices have samples in columns and genes/transcripts in rows

library(NOISeq)

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
  ## Files to run DE tool on ##
  subDirMatrices<-paste0("em/",expType)
  file.names <- dir(file.path(mainDir,subDirMatrices), pattern =".txt")
# loop through files #  
  for (fileName in file.names){
    RaEm<-unlist(strsplit(fileName, split='_', fixed=TRUE))[1]
    RaEmDe<-paste0(RaEm,"No")
    
    ## Define output directory and create if necessary
    subDirOut<-paste0(RaEmDe,"_",expType)
    if(file.exists(file.path(mainDirOut,subDirOut))){print(paste0("Combination ", RaEmDe,expType, " has been run already. Delete or rename ALL output folders for this comparison to rerun."))
    } else{
      print(paste0("Analyzing ", RaEmDe, "_", expType, "\n"))
      dir.create(file.path(mainDirOut,subDirOut))
    
    ## load in Expression Matrix ##
    expMat<-read.table(file.path(mainDir,subDirMatrices,fileName), header=TRUE, row.names=1)
    
    ## Set up Count Matrix and Groupings ##
    expMat[]<-as.integer(as.matrix(expMat))
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
    #### Run NOIseqBIO ####
    ######################

    ## create factors data.frame ##
    myfactors<-data.frame(group=groupings, row.names=colnames(expMat))
  
    ## Create a NOISeq object ##
    myData<-readData(data=expMat,factors=myfactors)
  
    ## Run NOISeqBIO, normalization only on counts data##, 
    if(expType=="countsGn"| expType=="countsTx"){noiSeqBioOut<-noiseqbio(myData, norm="tmm", factor="group", filter=1)}else{noiSeqBioOut<-noiseqbio(myData, norm="n", factor="group", filter=1)}

    ## Output Results
    results<-noiSeqBioOut@results[[1]]
    results$identifier<-row.names(results)
    names(results)[names(results)=="identifier"] <- identifier
    results<-merge(results,conversionTable,by= identifier)

    ## q is 1-FDR for NOIseqBIO
    noiSig = degenes(noiSeqBioOut, q = 0.95, M = NULL)


    ## Add annotations 
    if(identifier == "ENSG"){
      gnConvTable<-conversionTable[,c("ENSG","geneName")]
      noiSig<-merge(as.matrix(noiSig),gnConvTable, by.x="row.names", by.y="ENSG")
      noiSig<-unique(noiSig)
      colnames(noiSig) = c("ENSG",colnames(noiSig[2:length(noiSig)]))
    } else{
        noiSig<-merge(as.matrix(noiSig),conversionTable, by.x="row.names", by.y="ENST")
        colnames(noiSig) = c("ENST",colnames(noiSig[2:length(noiSig)]))
    }


    ######################
    #### Write output ####
    ######################
    
    if(compType == "gene"){
    write.table(unique(noiSig), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigGenes05.txt")),quote=FALSE)
    write.table(unique(results), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_allGenes.txt")),quote=FALSE)
    } else {
    write.table(noiSig, sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigTranscripts05.txt")),quote=FALSE)
    write.table(results, sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_allTranscripts.txt")),quote=FALSE)
    }
    write.table(unique(noiSig$geneName),sep="\t", row.names=FALSE,file.path(sigDirOut,expType,paste0(RaEmDe,"_","sigSymbol.txt")), quote=FALSE, col.names=FALSE)

    }
  }
}
