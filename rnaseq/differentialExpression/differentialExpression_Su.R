## Rscript to run Sleuth on all compatible expression outputs

library(sleuth)
library(stringr)

## Specify sample specific data

mainDir <- getwd()
controlGroup <- "classical"
testGroup <- "nonclassical"
alignCodes <- c("KaKa4Su", "SlSl4Su", "SfSn4Su", "SqSn4Su")
mainDirIn <- "mapAndModel"

mainDirOut <- "de/differentialExpression"
sigDirOut <- "de"

conversionTable <- read.table("conversion_tgs.txt", header=TRUE)

for (alignCode in alignCodes){
    RaEmDe = paste0(substring(alignCode,1,4), "Su")
    expType <- "tpmTx"

    ## Define output directory and create if necessary
    subDirOut<-paste0(RaEmDe,"_",expType)
    if(file.exists(file.path(mainDirOut,subDirOut))){print(paste0("Combination ", RaEmDe,expType, " has been run already. Delete or rename ALL output folders for this comparison to rerun."))
    } else{
      print(paste0("Analyzing ", RaEmDe, "_", expType, "\n"))
      dir.create(file.path(mainDirOut,subDirOut))
    
      ## Specify location of kallisto results and samples
      sample_id <- dir(file.path(mainDir, mainDirIn), pattern=paste0(controlGroup,"*|",testGroup,"*"))
      kal_dirs <- sapply(sample_id, function(id) file.path(mainDir, mainDirIn, id, alignCode))
      condition <- str_extract(sample_id, "[A-z]+")
      s2c <- data.frame(sample_id, condition)
      s2c <- dplyr::mutate(s2c, path=kal_dirs)
      colnames(s2c) = c("sample", "condition", "path")

      t2g = data.frame(conversionTable$ENSG, conversionTable$ENST, conversionTable$geneName)
      names(t2g) = c("ens_gene", "target_id", "geneName")

      ####################
      #### Run Sleuth ####
      ####################

      ## Create sleuth object and run differential expression testing
      so <- sleuth_prep(s2c, ~condition, target_mapping = t2g)
      so <- sleuth_fit(so)
      testCondition <- paste("condition", testGroup, sep="")
      so <- sleuth_wt(so, testCondition)

      ## Get results and identify significant transcripts at a q value of 0.05
      results_table <- sleuth_results(so, testCondition)
      sig <- subset(results_table, qval<0.05)

      ######################
      #### Write Output ####
      ######################

      comparison <- paste0(controlGroup,"VS",testGroup)
      write.table(results_table, sep="\t", row.names=FALSE, file.path(mainDirOut,subDirOut,paste0(RaEmDe,"_", expType, "_results.txt")), quote=FALSE)
      write.table(sig, sep="\t", row.names=FALSE, file.path(mainDirOut, subDirOut,paste0(RaEmDe, "_", expType, "_sigTranscripts05.txt")),quote=FALSE)
      write.table(unique(sig$geneName), sep="\t", row.names=FALSE, col.names=FALSE, file.path(sigGenesDir, paste0(RaEmDe,"_sigSymbol.txt")), quote=FALSE)
    }
}
