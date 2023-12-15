### Descriptive stats and statistical analysis
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
if(T){
  rm(list=ls())
  library(BayesFactor)
  library(ggplot2)
  library(stringr)
  library(lme4)
  source("functions/transformData.R")
  
  data=read.delim("dataset\\datasetOS.csv",sep=";")
  data_fp=read.delim("dataset\\datasetFP.csv",sep=";")
  data_rpe=read.delim("dataset\\datasetRPE.csv",sep=";")
  
  
  #outlier exclusion
  #get ids with performance below chance level
  below_chance_idData=getBelowChanceIdCond(data)
  # exclude participants with performance below chance level
  data = data[!data$ID %in% below_chance_idData,] 
  data_fp = data_fp[!data_fp$ID %in% below_chance_idData,]
  data_rpe = data_rpe[!data_rpe$ID %in% below_chance_idData,]
  
  #exclude rts more than 3 sds above or below mean
  data_fp = data_fp[!data_fp$outlier,]
  data = data[!data$outlier,]
  
  data$ID = as.factor(data$ID)
  data_fp$ID = as.factor(data_fp$ID)
  data_rpe$ID = as.factor(data_rpe$ID)
  
  data$deg = as.factor(data$deg)
  data_fp$deg = as.factor(data_fp$deg)
  
  data$block = as.factor(data$block)
  data_fp$block = as.factor(data_fp$block)
  data_rpe$block = as.factor(data_rpe$block)
}




########## Main Hypothesis ##########
##### Mean Amplitude #####
scaling_factor_fix = 0.5
set.seed(1234)
bf_full = anovaBF(MeanAmplitude ~ deg + block + deg:block + ID,
                  data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
rbind(row.names(extractBF(bf_full)), round(log10(extractBF(bf_full)$bf),3))





##### Sway Velocity #####
scaling_factor_fix = 0.5
set.seed(1234)
bf_full = anovaBF(SwayVelocity ~ deg + block + deg:block + ID,
               data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
rbind(row.names(extractBF(bf_full)), round(log10(extractBF(bf_full)$bf),3))



#### Max Range X ####
scaling_factor_fix = 0.5
set.seed(1234)
bf_full = anovaBF(MaxRangeX ~ deg + block + deg:block + ID,
                  data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
rbind(row.names(extractBF(bf_full)), round(log10(extractBF(bf_full)$bf),3))


#### Max Range Y ####
scaling_factor_fix = 0.5
set.seed(1234)
bf_full = anovaBF(MaxRangeY ~ deg + block + deg:block + ID,
                  data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
rbind(row.names(extractBF(bf_full)), round(log10(extractBF(bf_full)$bf),3))


########## Secondary Hypotheses ##########
##### Hypothesis S1: Cognitive Effort #####
data_cRPE = data_rpe[data_rpe$type=="cRPE",]

scaling_factor_fix = 0.5
set.seed(1234)
bf_block = anovaBF(value ~ block + ID, data = data_cRPE, whichRandom = "ID", rscaleFixed = scaling_factor_fix)

# t-tests of cRPE means between conditions
means_ment = c()
means_visu = c()
means_visu_stop = c()

for (ID in unique(data_cRPE$ID)) {
  means_ment = c(means_ment, mean(data_cRPE$value[data_cRPE$block=="mental" & data_cRPE$ID==ID]))
  means_visu = c(means_visu, mean(data_cRPE$value[data_cRPE$block=="visual" & data_cRPE$ID==ID]))
  means_visu_stop = c(means_visu_stop, mean(data_cRPE$value[data_cRPE$block=="visualStop" & data_cRPE$ID==ID]))
}

mean(means_ment)
mean(means_visu)
mean(means_visu_stop)

ttestBF(means_ment, means_visu, paired = TRUE)
ttestBF(means_ment, means_visu_stop, paired = TRUE)
ttestBF(means_visu, means_visu_stop, paired = TRUE)



##### Hypothesis S2: Physical Effort #####
data_RPE = data_rpe[data_rpe$type=="RPE",]

scaling_factor_fix = 0.5
set.seed(1234)
bf_block = anovaBF(value ~ block + ID, data = data_RPE, whichRandom = "ID", rscaleFixed = scaling_factor_fix)



##### Hypothesis S3: Reaction Time #####
data_rt = data[data$correct==1 & data$outlier==FALSE,]

scaling_factor_fix = 0.5
set.seed(1234)
bf_full = anovaBF(reactionTime ~ deg + block + deg:block + ID,
                  data = data_rt, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
rbind(row.names(extractBF(bf_full)), round(bf_full@bayesFactor$bf,3)) # bf_full@bayesFactor contains log10(BFs) without Inf problems

# t-tests of reaction time means between conditions
means_ment = c()
means_visu = c()
means_visu_stop = c()

for (ID in unique(data_rt$ID)) {
  means_ment = c(means_ment, mean(data_rt$reactionTime[data_rt$block=="mental" & data_rt$ID==ID]))
  means_visu = c(means_visu, mean(data_rt$reactionTime[data_rt$block=="visual" & data_rt$ID==ID]))
  means_visu_stop = c(means_visu_stop, mean(data_rt$reactionTime[data_rt$block=="visualStop" & data_rt$ID==ID]))
}

mean(means_ment)
mean(means_visu)
mean(means_visu_stop)

ttestBF(means_ment, means_visu, paired = TRUE)
ttestBF(means_ment, means_visu_stop, paired = TRUE)
ttestBF(means_visu, means_visu_stop, paired = TRUE)



##### Hypothesis S4: Pupil Size #####
pupildata = data[data$correct==1 & data$outlier==FALSE,]
pupildata = pupildata[!is.na(pupildata$PupilSizeDif),]

scaling_factor_fix = 0.5
set.seed(1234)
bf_full = anovaBF(PupilSizeDif ~ deg + block + deg:block + ID,
                  data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
rbind(row.names(extractBF(bf_full)), round(log10(extractBF(bf_full)$bf),3))

# t-tests of reaction time means between conditions
means_ment = c()
means_visu = c()
means_visu_stop = c()

for (ID in unique(pupildata$ID)) {
  means_ment = c(means_ment, mean(pupildata$PupilSizeDif[pupildata$block=="mental" & pupildata$ID==ID]))
  means_visu = c(means_visu, mean(pupildata$PupilSizeDif[pupildata$block=="visual" & pupildata$ID==ID]))
  means_visu_stop = c(means_visu_stop, mean(pupildata$PupilSizeDif[pupildata$block=="visualStop" & pupildata$ID==ID]))
}

mean(means_ment, na.rm=T)
mean(means_visu)
mean(means_visu_stop, na.rm=T)

ttestBF(means_ment[!is.na(means_ment) & !is.na(means_visu)], means_visu[!is.na(means_ment) & !is.na(means_visu)], paired = TRUE)
ttestBF(means_ment[!is.na(means_ment) & !is.na(means_visu_stop)], means_visu_stop[!is.na(means_ment) & !is.na(means_visu_stop)], paired = TRUE)
ttestBF(means_visu[!is.na(means_visu) & !is.na(means_visu_stop)], means_visu_stop[!is.na(means_visu) & !is.na(means_visu_stop)], paired = TRUE)
