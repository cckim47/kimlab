## Rscript to run EBseq on all counts expression matrices
## Expression matrices have samples in columns and genes/transcripts in rows

library(EBSeq)

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
        RaEmDe<-paste0(RaEm,"Eb")

        ## Define output directory and create if necessary
        subDirOut<-paste0(RaEmDe,"_",expType)
        if(file.exists(file.path(mainDirOut,subDirOut))){print(paste0("Combination ", RaEmDe,expType, " has been run already. Delete or rename ALL output folders for this comparison to rerun."))
        } else{
           print(paste0("Analyzing ", RaEmDe, "_", expType, "\n"))
           dir.create(file.path(mainDirOut,subDirOut))  

           ## Load in expression matrix 
           expMat<-as.matrix(read.table(file.path(mainDir,subDirMatrices,fileName), header=TRUE, row.names=1))
  
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
 
           ####################
           #### Run EBseq ####
           ###################

           groupString <- paste0(sum(groupings==group1),",",sum(groupings==group2))

           ## Generate library size vectors using median normalization
           Sizes=MedianNorm(expMat)

           ## Generate a vector of uncertainty groups for transcripts based on the number of transcripts per gene
           if (compType == "gene") {
              NgVector <- NULL
           }
           if (compType == "transcript") {
              NgList <- GetNg(conversionTable$ENST, conversionTable$ENSG)
              NgVector <- NgList$IsoformNgTrun
           }
  
           EBOut <- EBTest (Data = expMat, NgVector = NgVector, Conditions = groupings, sizeFactors = Sizes, maxround = 5)
           EBDERes <- GetDEResults(EBOut, FDR=0.05)
           resFC <- PostFC(EBOut)

           ## Merge gene name information into results
           if (compType == "gene"){
   	      resSymbol<-merge(as.matrix(EBDERes$DEfound),conversionTable[,c("ENSG","geneName")], by.x=1, by.y="ENSG")
              colnames(resSymbol) = c("ENSG",colnames(resSymbol[2:length(resSymbol)]))
  	   } else{
	      resSymbol<-merge(as.matrix(EBDERes$DEfound),conversionTable, by.x=1, by.y="ENST")
              colnames(resSymbol) = c("ENST",colnames(resSymbol[2:length(resSymbol)]))
	   }

           ######################
           #### Write output ####
           ######################

           write.table(cbind(EBDERes$PPMat,EBDERes$Status), sep="\t", file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_results.txt")), quote=FALSE, col.names=c(colnames(EBDERes$PPMat),"Status"), row.names=TRUE)
           write.table(cbind(resFC$PostFC, resFC$RealFC), sep="\t", file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_foldChanges.txt")), quote=FALSE, col.names=c("PostFC", "RealFC"), row.names=TRUE)
           write.table(EBDERes$DEfound, sep="\t", file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", comparison, "_DE.txt")), quote=FALSE, col.names=FALSE, row.names=FALSE)
           write.table(unique(resSymbol), sep="\t", file.path(mainDirOut, subDirOut, paste0(RaEmDe,"_",comparison, "_DE_Ann.txt")), quote=FALSE, col.names=TRUE, row.names=FALSE)
           write.table(unique(resSymbol$geneName), sep = "\t", file.path(sigDirOut, expType, paste0(RaEmDe,"_sigSymbol.txt")), quote=FALSE, col.names=FALSE, row.names=FALSE)
     
        }
    }
}
