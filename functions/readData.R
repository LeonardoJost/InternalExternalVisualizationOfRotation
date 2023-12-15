### Functions to read and modify data
#     Copyright (C) 2022  Leonardo Jost
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

source("functions/helpers.R")

library(plyr)
library(stringr)

#get questionnaireData
#verbose: detail of output
#folder: folder to search in for data
getQuestionnaireData=function(software,verbose,folder){
  if (verbose>1) {
    print("Reading questionnaire data from files ...")
  }
  if (software=="OSWeb"){
    questionnaireData=getDataOSWeb(verbose,folder,part="questionnaire")
    questionnaireData=subset(questionnaireData,select=c("ID","questionID","answer"))
    questionnaireData=reshape(questionnaireData, idvar = "ID", timevar = "questionID", direction = "wide")
    names(questionnaireData)=gsub("answer.","",names(questionnaireData))
    
  }else if (software=="OpenSesame"){
    questionnaireData=getQuestionnaireDataOpenSesame(verbose,folder, preText="", c("questionaire","questionnaire"),ending="csv")
  }else if (software=="jQuery"){
    questionnaireData=getQuestionnaireDataJQuery(verbose,folder, preText="",ending="csv")
  }
  if (verbose>1) {
    print(paste("Questionnaire data from",nrow(questionnaireData),"participants was read."))
  }
  return(questionnaireData)
}

#get questionnairedata for RPE questions (different block)
getRPEQuestionnaireData=function(software,verbose,folder){
  if (verbose>1) {
    print("Reading RPEquestionnaire data from files ...")
  }
  if (software=="OpenSesame"){
    questionnaireData=getRPEQuestionnaireDataOpenSesame(verbose,folder, preText="", c("mentalQuestionnaire","visualQuestionnaire","visualStopQuestionnaire"),ending="csv")
  }else {
    #not implemented
  }
  if (verbose>1) {
    print(paste("RPEQuestionnaire data from",nrow(questionnaireData),"participants was read."))
  }
  return(questionnaireData)
}

#get mental rotation data
#verbose: detail of output
#folder: folder to search in for data
#block: name of block of interest
getMRData=function(software,verbose,folder,block="main"){
  if (verbose>1) {
    print(paste("Reading mental rotation data for block",block,"from files"))
  }
  if (software=="OSWeb"){
    MRData=getDataOSWeb(verbose,folder,part=block)
  }else if (software=="OpenSesame"){
    MRData=getDataOpenSesame(verbose,folder,part=block)
  }
  if (verbose>1) {
    print(paste("Mental rotation data from",length(unique(MRData$ID)),"participants was read. (",nrow(MRData),"trials in total)"))
  }
  return(MRData)
}

#modifies the questionnairedata, calculates some additional information
#questionnaireData: dataset
modifyQuestionnaireData=function(questionnaireData,toFirstChars,toNums,cleanWhiteSpaces) {
  if (verbose>1) {
    print("Doing calculations on questionnaire data ...")
  }
  #transform values to numeric, remove white spaces, unify gender
  questionnaireData=cleanData(questionnaireData,toFirstChars,toNums,cleanWhiteSpaces)
  #rename columns to different names
  #colnames(questionnaireData) = make.unique(names(questionnaireData))
  if (verbose>1) {
    print("Calculations on questionnaire data finished.")
  }
  return(questionnaireData)
}

#modifies the mental rotation data, calculates some additional information
#verbose: detail of output
#MRData: dataset
#outlierFactor: trials deviating by more than outlierFactor*sd from mean will be classified as outliers
modifyMRData=function(verbose,MRData,outlierFactor) {
  if (verbose>1) {
    print("Doing calculations on mental rotation data ...")
  }
  #rename variables
  MRData$deg=toNumeric(MRData$angle)
  MRData$reactionTime=MRData$response_time
  MRData$angle=NULL
  MRData$response_time=NULL
  MRData$correctSide=MRData$correct_response
  MRData$modelNumber=paste("m",stringToNum(MRData$model),sep="")
  #save original degrees of rotation
  MRData$originalDegrees=MRData$deg
  #modify angles to 360-angle if angle>180, but keep information
  MRData$direction=ifelse(MRData$deg>180,"-",ifelse(MRData$deg==0 | MRData$deg==180,"0","+"))
  MRData$deg=ifelse(MRData$deg>180,360-MRData$deg,MRData$deg)
  MRData$meanPupilSize=MRData$pupilSizeSum/MRData$frames
  #mark outliers
  MRData=sortOutliers(verbose,MRData,outlierFactor)
  MRData$type=ifelse(MRData$correct==1,"hit","incorrect")
  MRData$typeOutlier=ifelse(MRData$outlier,paste(toChar(MRData$type),"Outlier",sep=""),toChar(MRData$type))
  #get maximum of all pupil sizes (pupilSize is last frame)
  MRData$maxPupilSize=do.call(pmax,list(MRData$maxPupilSize,MRData$pupilSize))
  #drop unnecessary rows
  MRData[,c("stimulusBottom","stimulusLeft","stimulusRight","correct_response","axis","model","pupilSize","zaehler_packetID","forceplate_startzeit")]=NULL
  return(MRData)
}

#get baseline pupil size and gaze data per trial
getBaselinePupilSizeAndGaze=function(verbose,showStimulusMRData) {
  #choose first 3 values for baseline pupil size
  showStimulusMRData$baselinePupilSizes=showStimulusMRData$pupilSize
  showStimulusMRData$baselinePupilSizes[showStimulusMRData$frames>3]=NA
  showStimulusMRData$baselinePupilSizes[showStimulusMRData$baselinePupilSizes<1]=NA
  #get focus on left, right, or alternative
  showStimulusMRData$gazeFocus=ifelse(showStimulusMRData$pupilSize<1,"",
                                      ifelse(showStimulusMRData$stimulusBottom==1,"t",
                                             ifelse(showStimulusMRData$stimulusLeft==1,"l",
                                                    ifelse(showStimulusMRData$stimulusRight==1,"r","x"))))  #x should not occur
  #rewrite left and right to correct and wrong alternative
  showStimulusMRData$gazeFocusCorrect=ifelse(showStimulusMRData$pupilSize<1,"",
                                             ifelse(showStimulusMRData$stimulusBottom==1,"t",
                                                    ifelse(showStimulusMRData$gazeFocus==substr(showStimulusMRData$correct_response,1,1),"c",
                                                           ifelse(showStimulusMRData$gazeFocus!="","w","x"))))
  #group data by stimulus
  showStimulusMRDataSummarized=ddply(showStimulusMRData,
                            .(ID,block,startTimeOfStimulus,startTimeOfBlock),
                            summarize,
                            #pupil size -> median of first 3
                            baselinePupilSize=median(baselinePupilSizes,na.rm=T),
                            #max pupil size
                            maxPupilSize=max(pupilSize,na.rm=T),
                            #order of fixations
                            gazeOrder=paste0(gazeFocus,collapse=""),
                            gazeOrderCW=paste0(gazeFocusCorrect,collapse=""))
  #number of measured fixations for trial
  showStimulusMRDataSummarized$lengthFocus=nchar(showStimulusMRDataSummarized$gazeOrder)
  #split gaze orders in halfs
  #showStimulusMRDataSummarized$gazeOrderFirstHalf=substr(showStimulusMRDataSummarized$gazeOrder,1,ceiling(showStimulusMRDataSummarized$lengthFocus/2))
  showStimulusMRDataSummarized$gazeOrderCWSecondHalf=substr(showStimulusMRDataSummarized$gazeOrderCW,floor(showStimulusMRDataSummarized$lengthFocus/2+1),showStimulusMRDataSummarized$lengthFocus)
  #remove length again
  showStimulusMRDataSummarized$lengthFocus=NULL
  #remove repetitions in gaze order
  showStimulusMRDataSummarized$gazeOrderSingles=unlist(gsub('([[:alpha:]])\\1+', '\\1', showStimulusMRDataSummarized$gazeOrder))
  #get unique orders
  showStimulusMRDataSummarized$gazeOrderUnique=unlist(lapply(strsplit(showStimulusMRDataSummarized$gazeOrder,""),function(x) paste0(unique(x),collapse = "")))
  #get proportion of correct and wrong alternative in second half
  showStimulusMRDataSummarized$gazeCWSecondHalfProportions=unlist(lapply(showStimulusMRDataSummarized$gazeOrderCWSecondHalf,function(x) str_count(x,"c")/(str_count(x,"w")+str_count(x,"c"))))
  #remove "ShowStimulus" for merging blocks later
  showStimulusMRDataSummarized$block=lapply(strsplit(showStimulusMRDataSummarized$block,split="Show"),`[[`,1)
  return(showStimulusMRDataSummarized)
}


#mark outliers, which deviate by more than sdFactor*sd from the mean by degree
sortOutliers=function(verbose,MRData,sdFactor) {
  degrees=levels(as.factor(MRData$deg))
  conditions=levels(as.factor(MRData$block))
  #allData$type=toChar(allData$type)
  MRData$outlier=FALSE
  for(degree in degrees) {
    for (condition in conditions) {
      degreeConditionSubset=MRData[which(MRData$deg==degree & MRData$block==condition),]
      meanRT=mean(degreeConditionSubset$reactionTime)
      sdRT=sd(degreeConditionSubset$reactionTime)
      if (verbose>2){
        print(paste("mean+-sd at angle ",degree," : ",meanRT,"+-",sdRT,sep=""))
      }
      MRData$outlier[which(MRData$deg==degree & MRData$block==condition & abs(MRData$reactionTime-meanRT)>sdFactor*sdRT)]=TRUE
    }
  }
  #allData$type=as.factor(allData$type)
  if (verbose>1) {
    print(paste(sum(MRData$outlier),"outliers detected (deviating by more than",
                outlierFactor,"standard deviations from mean (by degree)"))
  }
  return(MRData)
}


#reads data from files
#verbose: detail of output
#folder: folder to search in for files
#preText: Filter, only get files which start with preText
#part: Filter, only get parts of data in blocks whose name matches part
#ending: filetype of files
getDataOpenSesame=function(verbose, folder, preText="", part="main",ending="csv") {
  #get files in folger (Reaction Time Data)
  fileNames=getFileNames(folder,preText,ending)
  if (verbose>2) {
    print("list of files:")
    print(fileNames)
  }
  #initialize empty dataframe
  dat=data.frame()
  #loop through all files
  for (fileIndex in 1:length(fileNames)) {
    fileName=fileNames[fileIndex]
    #read data in file as table
    rawData=read.csv(paste(folder,fileName,sep=""),header=TRUE,fill=TRUE, sep=",")
    #choose only specified block
    dataset=rawData[rawData$aaBlock %in% part,]
    # dataset$numberInBlock=ave(dataset[,1],                 # Create numbering variable
    #                           dataset$aaBlock,
    #                           FUN = seq_along)
    if (verbose>3) {
      print(paste("read", nrow(dataset), "values from file:",fileName,"\n"))
    }
    #add to dataset
    dat=rbind(dat,dataset)
  }
  #change names
  dat$block=dat$aaBlock
  dat$ID=dat$aaID
  dat$aaBlock=NULL
  dat$aaID=NULL
  return(dat)
}

getRPEQuestionnaireDataOpenSesame=function(verbose, folder, preText="", part=c("mentalQuestionnaire","visualQuestionnaire","visualStopQuestionnaire"),ending="csv") {
  #get files in folger (Reaction Time Data)
  fileNames=getFileNames(folder,preText,ending)
  if (verbose>2) {
    print("list of files:")
    print(fileNames)
  }
  #initialize empty dataframe
  dat=data.frame()
  #loop through all files
  for (fileIndex in 1:length(fileNames)) {
    fileName=fileNames[fileIndex]
    #read data in file as table
    rawData=read.csv(paste(folder,fileName,sep=""),header=TRUE,fill=TRUE, sep=",")
    #choose only specified block
    dataset=rawData[rawData$aaBlock %in% part,]
    if (verbose>3) {
      print(paste("read", nrow(dataset), "values from file:",fileName,"\n"))
    }
    #add to dataset
    dat=rbind(dat,dataset)
  }
  #change names and keep only 4 columns
  names(dat)[1:4]=c("block","ID","value","type")
  dat[,5:ncol(dat)]=NULL
  return(dat)
}

#reads data from OSWeb/JATOS files
#verbose: detail of output
#folder: folder to search in for files
#preText: Filter, only get files which start with preText
#part: Filter, only get parts of data in blocks whose name matches part
#ending: filetype of files
getDataOSWeb=function(verbose, folder, preText="", part="main",ending="csv") {
  #get files in folger (Reaction Time Data)
  fileNames=getFileNames(folder,preText,ending)
  if (verbose>2) {
    print("list of files:")
    print(fileNames)
  }
  #initialize empty dataframe
  dat=data.frame()
  #loop through all files
  for (fileIndex in 1:length(fileNames)) {
    fileName=fileNames[fileIndex]
    #read data in file as table
    rawData=read.csv(paste(folder,fileName,sep=""),header=TRUE,fill=TRUE, sep=",")
    #choose only specified block
    dataset=rawData[rawData$aaBlock %in% part,]
    if (verbose>3) {
      print(paste("read", nrow(dataset), "values from file:",fileName,"\n"))
    }
    #add to dataset
    dat=rbind(dat,dataset)
  }
  #change names
  dat$block=dat$aaBlock
  dat$ID=dat$aaID #aaID for old version
  dat$aaBlock=NULL
  dat$aaID=NULL
  return(dat)
}

#modified 
getQuestionnaireDataJQuery=function(verbose, folder, preText="",ending="csv") {
  #get files in folger (Reaction Time Data)
  fileNames=getFileNames(folder,preText,ending)
  if (verbose>2) {
    print("list of files:")
    print(fileNames)
  }
  #initialize empty dataframe
  dat=data.frame()
  #loop through all files
  for (fileIndex in 1:length(fileNames)) {
    fileName=fileNames[fileIndex]
    #read data in file as table
    rawData=read.csv(paste(folder,fileName,sep=""),header=TRUE,fill=TRUE, sep=",")
    #choose only specified block
    dataset=rawData[rawData$aaBlock=="",]
    if (verbose>3) {
      print(paste("read", nrow(dataset), "values from file:",fileName,"\n"))
    }
    #add to dataset
    dat=rbind(dat,dataset)
  }
  #change names
  #dat$block="questionnaire"
  dat$ID=dat$workerId #aaId for old version
  dat$aaBlock=NULL
  dat$workerId=NULL
  return(dat)
}

#reads data from files
#verbose: detail of output
#folder: folder to search in for files
#preText: Filter, only get files which start with preText
#part: Filter, only get parts of data in blocks whose name matches part
#ending: filetype of files
getQuestionnaireDataOpenSesame=function(verbose, folder, preText="", part=c("questionaire","questionnaire"),ending="csv") {
  #get files in folder
  fileNames=getFileNames(folder,preText,ending)
  if (verbose>2) {
    print("list of files:\n")
    print(fileNames)
  }
  #initialize empty dataframe
  dat=data.frame()
  #loop through all files
  for (fileIndex in 1:length(fileNames)) {
    fileName=fileNames[fileIndex]
    #read data in file as table
    rawData=read.csv(paste(folder,fileName,sep=""),header=TRUE,fill=TRUE, sep=",")
    #choose only specified block
    dataset=rawData[rawData$aaBlock %in% part,]
    #add interesting data to vector 
    values=append(toChar(dataset[,3]),dataset$aaID[1])
    if (verbose>3) {
      print(paste("read values for file:",fileName,"\n"))
      print(values)
    }
    #add to dataset
    dat=rbind(dat,values,stringsAsFactors = FALSE)
    #set names according to questionIDs
    if (fileIndex==1) {
      names(dat)=append(toChar(dataset[,4]),"ID")
    }
  }
  return(dat)
}