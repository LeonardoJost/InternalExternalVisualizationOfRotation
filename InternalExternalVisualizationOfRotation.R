### main program
#     Copyright (C) 2023  Leonardo Jost and Markus Siebertz
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
source("functions/readData.R", encoding="utf-8")
source("functions/generateGraphsAndTables.R", encoding="utf-8")
source("functions/transformData.R")

##create output directories, if they don't exist (outputs warnings otherwise)
dir.create("figs")
dir.create("output")

##options, parameters
options(digits=6)
software=("OpenSesame") #OSWeb or OpenSesame (classic)
#set data folder
folder="Logfiles\\"
verbose=3 #detail of output
questionnaireOutFile="output\\questionnaire" #.csv added at end, leave empty if no output desired
outlierFactor=3 #factor of sd to define outliers in MR
block=c("mental","visual","visualStop")#name of interesting block of data
showStimulusBlocks=c("mentalShowStimulus","visualShowStimulus","visualStopShowStimulus")
questionnaireDataColsDescriptive=c("ID","Gender","Age","Experience") #which questionnaire columns shall be kept for descriptive reports
questionnaireDataCols=c("ID") #which questionnaire columns shall be kept for statistical analysis

##read and write data
#read data
questionnaireData=getQuestionnaireData(software,verbose,folder)
questionnaireDataRPE=getRPEQuestionnaireData(software,verbose,folder)
MRData=getMRData(software,verbose,folder,block)
showStimulusMRData=getMRData(software,verbose,folder,showStimulusBlocks)
#merge summarized data from ShowStimulus with MRData
#calculate baseline pupil size first
showStimulusMRDataSummarized=getBaselinePupilSizeAndGaze(verbose,showStimulusMRData)
MRData=merge(MRData,showStimulusMRDataSummarized,by=c("ID","startTimeOfStimulus","startTimeOfBlock","block"))
#modify data 
questionnaireData=modifyQuestionnaireData(questionnaireData,c("Gender"),c("Age"),c())
MRData=modifyMRData(verbose,MRData,outlierFactor) 
#calculate means from questionnaire (and save to csv)
calculateMeansQuestionnaire(verbose,questionnaireData,questionnaireOutFile,"")
#remove not relevant questionnaire data before merging
questionnaireData=subset(questionnaireData,select=questionnaireDataColsDescriptive)
#unify data
datasetDesc=merge(MRData,questionnaireData,by="ID")

#remove not analyzed questionnaire data to protect participant identity
questionnaireData=subset(questionnaireData,select=questionnaireDataCols)
#unify data
dataset=merge(MRData,questionnaireData,by="ID")
#read forceplate data
datasetForceplate=read.delim("Logfiles\\InternalExternalVisualization_ForcePlateData157_trialsdeleted.txt", sep = ",")
#merge with force plate data
data_list = getData(dataset,datasetForceplate,questionnaireDataRPE,157)
#get subsets of data_list
#for analysis
data = data_list[[1]]
data_fp = data_list[[2]]
data_rpe = data_list[[3]]

#anonymise IDs to protect participant identity
#convert to factors and use same levels for all
data$ID=factor(data$ID)
data_fp$ID=factor(data_fp$ID,levels=levels(data$ID))
data_rpe$ID=factor(data_rpe$ID,levels=levels(data$ID))

#random permutation
randPerm=sample.int(length(levels(data$ID)))
#randomize all levels with the same permutation
levels(data$ID)=paste("id",randPerm,sep="")
levels(data_fp$ID)=paste("id",randPerm,sep="")
levels(data_rpe$ID)=paste("id",randPerm,sep="")

#randomize order of data frames
data=data[sample(1:nrow(data)),]
data_fp=data_fp[sample(1:nrow(data_fp)),]
data_rpe=data_rpe[sample(1:nrow(data_rpe)),]


#save full randomized dataset to csv
write.table(data,file="output\\datasetOS.csv",sep=";", row.names = F)
write.table(data_rpe,file="output\\datasetRPE.csv",sep=";", row.names = F)
write.table(data_fp,file="output\\datasetFP.csv",sep=";", row.names = F)


