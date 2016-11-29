## Program to convert bootstrap estimates from salmon and sailfish into output ready for sleuth
## Uses program wasabi
## In command line call, supply RA/EM combination, sample base name, average fragment length, average fragment length SD

library(wasabi)

## Read command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4){
        stop("Must supply subdirectory as a command line argument.  Ex: SlSl4Su \n
		Must also supply sample base name.  Ex: classical
		Also supply average fragment length.  Ex: 360
		Also supply fragment length sd. Ex: 200", call.=FALSE)
		
        }

subdir <- args[1]
sampleName <- args[2]
fragLength <- args[3]
fragSD <- args[4]

## Find samples
sample_id <- dir(pattern=paste0("^",sampleName,"*"))
sfdirs <- file.path(c(sample_id),subdir)

## Prepare for sleuth
prepare_fish_for_sleuth(sfdirs, fallback_mu=fragLength, fallback_sd=fragSD)
