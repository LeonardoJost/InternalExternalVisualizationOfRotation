### Outlier exclusion
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



#script for possible outlier exclusion -> lead to no exclusion

# ##### Check for missing gaze points for stim and/or both targets #####
#check performance in trials with no gaze points on stim or both targets
sum(data$stimulusBottomSum==0)
sum(data$stimulusLeftSum==0 & data$stimulusRightSum==0)
#which(data$stimulusBottomSum==0)
#which(data$stimulusLeftSum==0 & data$stimulusRightSum==0)
mean(data$correct)
mean(data$correct[data$stimulusBottomSum==0])
mean(data$correct[data$stimulusLeftSum==0 & data$stimulusRightSum==0])

# check performance in trials with no gaze points on stim or both targets where eye tracker worked
data_pup_no_na = data[!is.na(data$PupilSizeDif),]
sum(data_pup_no_na$stimulusBottomSum==0)
sum(data_pup_no_na$stimulusLeftSum==0 & data_pup_no_na$stimulusRightSum==0)
sum((data_pup_no_na$stimulusLeftSum==0 & data_pup_no_na$stimulusRightSum==0) | data_pup_no_na$stimulusBottomSum==0)
#which(data_pup_no_na$stimulusBottomSum==0)
#which(data_pup_no_na$stimulusLeftSum==0 & data_pup_no_na$stimulusRightSum==0)
mean(data_pup_no_na$correct)
mean(data_pup_no_na$correct[data_pup_no_na$stimulusBottomSum==0])
mean(data_pup_no_na$correct[data_pup_no_na$stimulusLeftSum==0 & data_pup_no_na$stimulusRightSum==0])

#aggregate by ids
dataTMP=data_pup_no_na
dataTMP$noLookAlternative=ifelse(dataTMP$stimulusLeftSum==0 & dataTMP$stimulusRightSum==0,T,F)
dataTMP$noLookTarget=ifelse(dataTMP$stimulusBottomSum==0,T,F)
noLooksById=ddply(dataTMP,.(ID),summarize,
                  noLookAlts=sum(noLookAlternative),
                  noLookTgts=sum(noLookTarget),
                  meanAcc=mean(correct),
                  meanAccNoAlt=weighted.mean(correct,noLookAlternative),
                  meanAccNoTgt=weighted.mean(correct,noLookTarget))
#do not exclude trials because average accuracy is very high indicated that the figures were observed


##### check for missing degrees and conditions #####
missing_values = function(data){
  missing_degree = FALSE
  missing_IDs_deg = c()
  missing_cond = FALSE
  missing_IDs_cond = c()
  missing_cond_deg = FALSE
  missing_IDs_cond_deg = c()
  for (ID in unique(data$ID)) {
    print(ID)
    print("trials per condition-degree combination: ")
    print(table(data$block[data$ID==ID],data$deg[data$ID==ID]))
    tmp_table_deg = table(data$deg[data$ID==ID])
    if(length(tmp_table_deg) < 4){
      missing_degree = TRUE
      
      print(paste("!!! Missing data for at least one whole degree for ID", ID))
      #print(tmp_table_deg)
      missing_IDs_deg = c(missing_IDs_deg, ID)
    }
    tmp_table_cond = table(data$block[data$ID==ID])
    if(length(tmp_table_cond) < 3){
      missing_cond = TRUE
      print(paste("!!! Missing data for at least one whole condition for ID", ID))
      #print(tmp_table_cond)
      missing_IDs_cond = c(missing_IDs_cond, ID)
    }
    tmp_table_cond_deg = table(data$block[data$ID==ID],data$deg[data$ID==ID])
    if(any(tmp_table_cond_deg==0)){
      print(paste("!!! Missing data for at least one specific condition-degree combination for ID", ID))
      missing_IDs_cond_deg = c(missing_IDs_cond_deg, ID)
      missing_cond_deg = TRUE
    }
    print("-----------------------------")
  }
  if (missing_degree){
    print("!!! Missing data for at least one whole degree for the following IDs:")
    print(missing_IDs_deg)
  }else{
    print("No missing data for any whole degree for any ID")
  }
  if (missing_cond){
    print("!!! Missing data for at least one whole condition for the following IDs:")
    print(missing_IDs_cond)
  }else{
    print("No missing data for any whole conditions for any ID")
  }
  if(missing_cond_deg){
    print("!!! Missing data for at least one specific condition-degree combination for the following IDs:")
    print(missing_IDs_cond_deg)
  }else{
    print("No missing data for any specific condition-degree combination for any ID")
  }
}
missing_values(data) # OpenSesame data
missing_values(data_fp) # force plate data
missing_values(data[data$correct==1,])
#only one missing data point (after chance level exclusion) -> impute

