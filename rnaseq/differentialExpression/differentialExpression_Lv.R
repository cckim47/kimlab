## Rscript to run limma + VST on all counts expression matrices
## Expression matrices have samples in columns and genes/transcripts in rows

library(limma)
library(plyr)
library(edgeR)
library(DESeq2)

########################
#### Wrapper  Input ####
########################

conversionTable<-read.table("conversion_tgs.txt", header=TRUE)

# define input and output directories
mainDir<-getwd()
mainDirOut<-"de/differentialExpression"
sigDirOut<-"de"

# define which expression types to run DE on #
expTypes <- c("countsGn","countsTx")

for (expType in expTypes){
  
  subDirMatrices<-paste0("em/",expType)
  file.names <- dir(file.path(mainDir,subDirMatrices), pattern =".txt")
  
  for (fileName in file.names){
    RaEm<-unlist(strsplit(fileName, split='_', fixed=TRUE))[1]
    RaEmDe<-paste0(RaEm,"Lv")
    
    ## Define output directory and create if necessary
    subDirOut<-paste0(RaEmDe,"_",expType)
    if(file.exists(file.path(mainDirOut,subDirOut))){print(paste0("Combination ", RaEmDe,expType, " has been run already. Delete or rename ALL output folders for this comparison to rerun."))
    } else{
      print(paste0("Analyzing ", RaEmDe, "_", expType, "\n"))
      dir.create(file.path(mainDirOut,subDirOut))
      
      ## Load in expression matrix 
      expMat<-read.table(file.path(mainDir,subDirMatrices,fileName), header=TRUE, row.names=1)
      
      ## Set up expression matrix and groupings 
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
      #### Run Limma + VST ####
      ######################
       
       
      ### Filter out lowly expressed genes - genes that are expressed in <= 17 samples (expressed in half the samples)
      cutoffCounts<-0
      expMat$sum<-rowSums(expMat>cutoffCounts)
      ## Define cutoff: 1/2 the number of samples ##
      cutoff<-(ncol(expMat)-1)/2
      ## remove genes that are not expressed in at least half the samples ##
      expMat<-subset(expMat, sum>=cutoff)
      ## remove sums column ##
      expMat<-expMat[-length(expMat)]

      ### Set up for limma VST ###       
      col3<-colnames(expMat)
      col3<-substr(col3,1,nchar(col3)-2)
      designMat<-data.frame(col1=1,col2=0, col3=col3,row.names=as.character(colnames(expMat)))
      designMat<-rename(designMat,c(col1="reference",col2=paste0(group1,"v",group2)))
      designMat[designMat$col3 == group1, paste0(group1,"v",group2)] = 1
      designMat<-designMat[-3]
       
      ## Limma VST ##
       
      colData<-data.frame(sample=as.character(colnames(expMat)), group=col3)
      row.names(colData)<-colData$sample
      expMatvst<-as.matrix(expMat)
      expMatvst[]<-as.integer(expMatvst)
      dds<-DESeqDataSetFromMatrix(countData=expMatvst, colData=colData, design= ~group)
      vsd<-varianceStabilizingTransformation(dds, blind=FALSE)
       
      fitVST<-lmFit(assay(vsd),designMat)
      fitVST<-eBayes(fitVST)
      resVST<-topTable(fitVST, coef=comparison,adjust="BH", n=Inf)
      resVST$UNIQID<-row.names(resVST)

      ## define if run on gene or transcript and merge with gene symbol ##
      names(resVST)[names(resVST)=="UNIQID"] <- identifier
      resSymbol<-merge(resVST,conversionTable, by=identifier)
      sig<-subset(resSymbol,adj.P.Val<0.05)

      ### Create Output Files ###
      if(compType == "gene"){
	      write.table(unique(sig), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigGenes05.txt")),quote=FALSE)
	      write.table(unique(resSymbol), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_allGenes.txt")),quote=FALSE)
	    } else {
	      write.table(sig, sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigTranscripts05.txt")),quote=FALSE)
        write.table(resSymbol, sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_allTranscripts.txt")),quote=FALSE)
        	}
      write.table(unique(sig$geneName),sep="\t", row.names=FALSE,file.path(sigDirOut,expType,paste0(RaEmDe,"_","sigSymbol.txt")), quote=FALSE, col.names=FALSE)

    }
  }
}

