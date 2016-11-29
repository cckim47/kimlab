## Rscript to run NBPSeq on all counts expression matrices
## Expression matrices have samples in columns and genes/transcripts in rows

library(NBPSeq)

########################
#### Wrapper  Input ####
########################

## Set up input and output files and directories

mainDir<-getwd()
mainDirOut<-"de/differentialExpression"
sigDirOut <- "de"

conversionTable<-read.table("conversion_tgs.txt", header=TRUE)

expTypes <- c("countsTx", "countsGn")

for (expType in expTypes){

    subDirMatrices<-paste0("em/",expType)
    file.names <- dir(file.path(mainDir,subDirMatrices), pattern =".txt")

    for (fileName in file.names){
        RaEm<-unlist(strsplit(fileName, split='_', fixed=TRUE))[1]
        RaEmDe<-paste0(RaEm,"Nb")
        
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
  
          #################### 
          #### Run NBPSeq ####
          ####################

          norm.factors <- estimate.norm.factors(expMat)
          res <- nbp.test(round(as.matrix(expMat)), groupings, group1, group2, norm.factors=norm.factors, print.level=5)
          resTable <- cbind(res$expression.levels, res$p.values, res$q.values)
          colnames(resTable) <- c(group1, group2, "pooled", "pval","qval")

          ## Merge gene name information into results
          if(identifier == "ENSG"){
	     resSymbol<-merge(resTable,conversionTable[,c("ENSG","geneName")], by.x="row.names", by.y=identifier)
	  } else{
	     resSymbol<-merge(resTable,conversionTable, by.x="row.names", by.y=identifier)
	  }

          ## Find significant genes or transcripts at q value of 0.05, removing those for which there was no test
          sig <- resSymbol[c(resSymbol$qval<0.05),]
          sig <- sig[c(!is.na(sig$qval)),]

          ######################
          #### Write output ####
          ######################

         if(identifier == "ENSG"){
	    write.table(unique(sig), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigGenes05.txt")),quote=FALSE)
	    write.table(unique(resSymbol), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_allGenes.txt")),quote=FALSE)
	 } else{
	    write.table(unique(sig), sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_sigTranscripts05.txt")),quote=FALSE)
            write.table(resSymbol, sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_allTranscripts.txt")),quote=FALSE)
	 }
         write.table(unique(sig$geneName),sep="\t", row.names=FALSE,file.path(sigDirOut,expType,paste0(RaEmDe,"_","sigSymbol.txt")), quote=FALSE, col.names=FALSE)

      }
    }
}
