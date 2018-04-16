#### Randomly subset fastq files to specified read depths using the ShortRead package. Modify desired read depth(s) in code below.

###################
#                 #
# Load Packages   #
#                 #
###################

library(ShortRead)

################################
#                              #
# set read depths (modify this)#
#                              #
################################


desiredReadDepths<-c(2e+10, 5e+10)


##########################################
#                                         #
#                                         #
#      Create Desired Folder Structure    #
#                                         #
#                                         #
###########################################

### 1, Create directories to place newly generated fastq files at each desired read depth (named by depth).
### 2. Identify all .fastq files within the working directory.  
### 3. Within each read depth directory, check for directory that shares name of the .fastq file to be sampled - create new subdirectory. (If subdirectory exists, will move on to next .fastq file to process)
### 4. Check that original .fastq file has enough reads to be sampled at desired read depth
### 5. Randomly sample reads from .fastq file - newly generated file will be named with the following convention /deiredReadDepth/originalFastqName/originalFastqName.fastq

mainDir<-getwd()
fastqRaw<- dir(pattern =".fastq")

for(i in desiredReadDepths){
  ## Run through desired sampling depths, check if these have been previously run ##
  if(file.exists(file.path(mainDir,i))){print(paste("FYI, the folder",i,"already exists."))} 
  else{
    ### Create a directory for output with i reads###
    print(paste0("Generating fastq files with",i, "reads"))
    dir.create(file.path(mainDir,i))
  }

    for (j in fastqRaw){
      ### keep sample name from original .fastq file
      sampleName<-unlist(strsplit(j,".",fixed=TRUE))[1]
      
      ### Check if output directory: /desiredReadDepth/sampleName exists - if it exists, skip####
      outputDir<-file.path(mainDir,i,sampleName)
      if(file.exists(outputDir)){print(paste("The fastq file for",sampleName,"in",n,"already exists. Delete or rename directory to regenerate this sample at this depth."))} 
      else{
        print(paste0("Processing ", sampleName, "with ", n, "reads"))
        ### Make the output directory with the sample name ###
        dir.create(outputDir)
      
      ###################################################
      #                                                 #
      #  create subsetted fastq file using shortread    #
      #                                                 #
      ###################################################
      sampled<-FastqSampler(j, i)
      extract<-yield(sampled)
      
      ### Check that there were enough reads in the the original file ###
      if(length(extract) < i){
        print(paste0("not enough reads in orginal ", j," file"))}
      else{
        writeFastq(extract, file.path(outputDir,paste0(sampleName,".fastq")), mode="w", compress= FALSE, full=TRUE)
      }
    }
  }
  closeAllConnections()
}