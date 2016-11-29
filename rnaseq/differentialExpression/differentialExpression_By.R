## Rscript to run ballgown on all counts expression matrices
## Expression matrices have samples in columns and genes/transcripts in rows

library(baySeq)
library(parallel)

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
    RaEmDe<-paste0(RaEm,"By")
    
    ## Define output directory and create if necessary
    subDirOut<-paste0(RaEmDe,"_",expType)
    if(file.exists(file.path(mainDirOut,subDirOut))){print(paste0("Combination ", RaEmDe,expType, " has been run already. Delete or rename ALL output folders for this comparison to rerun."))
    } else{
      print(paste0("Analyzing ", RaEmDe, "_", expType, "\n"))
      dir.create(file.path(mainDirOut,subDirOut))
      print(paste("directory created for", RaEmDe))
      
      
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
      
      ##########################################
      #### Filter out Lowly Expressed Genes ####
      ##########################################
    
      ### Filter out lowly expressed genes - genes that are expressed in <= 17 samples (expressed in half the samples)
      cutoffCounts<-0
      expMat$sum<-rowSums(expMat>cutoffCounts)
      ## Define cutoff: 1/2 the number of samples ##
      cutoff<-(ncol(expMat)-1)/2
      ## remove genes that are not expressed in at least half the samples ##
      expMat<-subset(expMat, sum>=cutoff)
      ## remove sums column ##
      expMat<-expMat[-length(expMat)]

      
      ######################
      #### Run BaySeq ####
      ######################
         
      ## baySeq Setup         
      ### Detect available cores ###
      availCores<-detectCores(logical=FALSE)
      availCores<-availCores-2
      cl<-makeCluster(availCores)
         
      ## create lists for NDE and DE ##
      expMat<-as.matrix(expMat)
      NDE<-factor(matrix(1, nrow=1,ncol=length(groupings)))
      DE<-groupings
      groups<-list(NDE=NDE,DE=DE)
         
      ## Create a baySeq object ##
      CD<-new("countData",data=expMat, replicates=groupings,groups=groups)
      ##infer library sizes - also the option of supplying this ##
      libsizes(CD)<-getLibsizes(CD, estimationType="quantile")
         
      ## Neg Binomial, run on suggested 10000 iterations - for testing recommend reducing to 1000 ##
      CD<-getPriors.NB(CD, samplesize=10000, estimation = "QL",cl=cl)
	    CD<-getLikelihoods(CD, cl=cl)
      ### Save CD file to output directory ##
	    save(CD, file=file.path(mainDirOut,subDirOut,"CDposterior.RData"))
   	  results<-topCounts(CD, group="DE", number=length(expMat))
      results$identifier<-row.names(results)
         
      ### Merge Gene Names to Results ###
      names(results)[names(results)=="identifier"] <- identifier
      results<-merge(results,conversionTable,by= identifier)
      ## Define significance cutoff and subset sig genes ##
      sig<-subset(results, FDR.DE <0.05)
         
         
      ######################
      #### Write output ####
      ######################
      
         if(compType == "gene"){
           write.table(unique(sig), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigGenes05.txt")),quote=FALSE)
           write.table(unique(results), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_allGenes.txt")),quote=FALSE)
         }
         else {
           write.table(sig, sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigTranscripts05.txt")),quote=FALSE)
           write.table(results, sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_allTranscripts.txt")),quote=FALSE)
         }
         write.table(unique(sig$geneName),sep="\t", row.names=FALSE,file.path(sigDirOut,expType,paste0(RaEmDe,"_","sigSymbol.txt")), quote=FALSE, col.names=FALSE)
         
    }
  }
}
