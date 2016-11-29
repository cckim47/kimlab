## Rscript to run BitSeq on all HsBs and KaBs expression estimates
## Does not take expression matrices

library(BitSeq)

#######################
#### Wrapper Input ####
#######################

## Set up input and output files and directories

mainDir <- getwd()
mainDirIn <- "mapAndModel/"
alignCodes <- c("HsBs", "KaBs")
controlGroup <- "classical"
testGroup  <- "nonclassical"

mainDirOut <- "de/differentialExpression"
sigGenesDir <- "de/fpkmTx"

conversionTable <- read.table("conversion_tgs.txt", header=TRUE)

for (RaEm in alignCodes) {

     RaEmDe <- paste0(RaEm,"Bs")

     ## Define output directory and create if necessary
     subDirOut<-paste0(RaEmDe,"_",expType)
     if(file.exists(file.path(mainDirOut,subDirOut))){print(paste0("Combination ", RaEmDe,expType, " has been run already. Delete or rename ALL output folders for this comparison to rerun."))
     } else{
        print(paste0("Analyzing ", RaEmDe, "_", expType, "\n"))
        dir.create(file.path(mainDirOut,subDirOut))

        ## Specify location of expression files
        controlID <- dir(file.path(mainDir,mainDirIn),pattern=paste0("^",controlGroup,"*"))
        controlFiles <- file.path(mainDir,mainDirIn, controlID, RaEm, "alignedExpressionRPKM.rpkm")
        testID <- dir(file.path(mainDir, mainDirIn),pattern=paste0("^",testGroup,"*"))
        testFiles <- file.path(mainDir, mainDirIn, testID,RaEm,"alignedExpressionRPKM.rpkm")

        ####################
        #### Run BitSeq ####
        ####################

        de1 <- getDE(list(controlFiles,testFiles), samples=FALSE, outPrefix=paste0(mainDirOut,"/", subDirOut, "/bitseqOut"), trInfoFile=file.path(controlID[1],RaEm,"alignedExpressionRPKM.tr"))

        ## Identify significant transcripts at PPLR (probability of positive log ratio) over 0.975 or below 0.025
        sigTranscripts = de1$pplr[abs(0.5-de1$pplr$pplr)>0.475,]
        colnames(sigTranscripts) = c("pplr","log2FoldChange","LowerCI","UpperCI","logMeanExpressionC1","logMeanExpressionC2")

        ## Merge gene name information into results
        sigGenes = merge(sigTranscripts,conversionTable, by.x="row.names", by.y="ENST", all.x=TRUE)[,c("Row.names","pplr","ENSG","geneName")]

        ######################
        #### Write Output ####
        ######################

        write.table(sigTranscripts, sep="\t",row.names=TRUE, file=paste0(mainDirOut,"/", subDirOut, "/sigTranscripts.txt"), quote=FALSE)
        write.table(sigGenes, sep="\t",row.names=FALSE, file=paste0(mainDirOut,"/", subDirOut, "/sigGenes.txt"), quote=FALSE)
        write.table(unique(sigGenes[,4]),row.names=FALSE, col.names=FALSE, file=paste0(sigGenesDir,"/",RaEmDe, "_sigSymbol.txt"), quote=FALSE) 
     }
}
