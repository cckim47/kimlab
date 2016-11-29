## Rscript to run DESeq2 on all count expression matrices
## Expression matrices have samples in columns and genes/transcripts in rows
library(DESeq2)

########################
#### Wrapper  Input ####
########################

## Set up input and output files and directories

mainDir<-getwd()
mainDirOut<-"de/differentialExpression"
sigDirOut<-"de"

conversionTable<-read.table("conversion_tgs.txt", header=TRUE)

expTypes <- c("countsGn","countsTx")


for (expType in expTypes){
  
  subDirMatrices<-paste0("em/",expType)
  file.names <- dir(file.path(mainDir,subDirMatrices), pattern =".txt")
  
  for (fileName in file.names){
    RaEm<-unlist(strsplit(fileName, split='_', fixed=TRUE))[1]
    RaEmDe<-paste0(RaEm,"De")
    
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
      #### Run DESeq2   ####
      ######################
      rounded <- round(expMat, 0)
      frame <- data.frame(groupings, row.names = colnames(expMat))
      dds<-DESeqDataSetFromMatrix(countData = rounded, colData = frame, design = ~ groupings)
      dds<-DESeq(dds)
    
      ## Generate results - alpha set to 0.05 to optimize filtering 
      res<-results(dds, alpha=0.05)
      resOrdered<-res[order(res$padj),]

      ######################
      #### Write output ####
      ######################
      if(identifier == "ENSG"){
	      resSymbol<-merge(as.matrix(resOrdered),conversionTable[,c("ENSG","geneName")], by.x="row.names", by.y="ENSG")
        colnames(resSymbol) = c("ENSG",colnames(resSymbol[2:length(resSymbol)]))
	      } else{
	       resSymbol<-merge(as.matrix(resOrdered),conversionTable, by.x="row.names", by.y="ENST")
        colnames(resSymbol) = c("ENST",colnames(resSymbol[2:length(resSymbol)]))
	      }
      sig <- resSymbol[c(resSymbol$padj<0.05),]
      sig <- sig[c(!is.na(sig$padj)),]
      write.table(resOrdered, sep="\t", file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_results.txt")), quote=FALSE, col.names=TRUE, row.names=TRUE)
      write.table(unique(resSymbol), sep="\t", file.path(mainDirOut, subDirOut, paste0(RaEmDe,"_",comparison, "_resultsAnn.txt")), quote=FALSE, col.names=TRUE, row.names=FALSE)
      write.table(unique(sig), sep = "\t", file.path(mainDirOut, subDirOut, paste0(RaEmDe, "_", comparison, "_adjp0.05.txt")), quote=FALSE, col.names=TRUE, row.names=FALSE)
      write.table(unique(sig$geneName), sep = "\t", file.path(sigDirOut, expType, paste0(RaEmDe,"_sigSymbol.txt")), quote=FALSE, col.names=FALSE, row.names=FALSE)
    }
  }
}
 
