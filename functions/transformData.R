### Merging Mental Rotation and Force Plate Data
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

load_and_prep_data = function(data,data_fp,data_rpe){
  
  ##### remove outlier based on rotation reaction time in opensesame data ##### 
  #data = data[!data$outlier,] 
  
  data$PupilSizeDif = data$meanPupilSize - data$baselinePupilSize
  data$maxPupilSizeDif=data$maxPupilSize-data$baselinePupilSize
  
  
  names(data)[names(data)=="aaBlock"] = "block"
  names(data)[names(data)=="aaID"] = "ID"
  
  names(data_fp)[names(data_fp)=="aaBlock"] = "block"
  names(data_fp)[names(data_fp)=="aaID"] = "ID"
  names(data_fp)[names(data_fp)=="angle"] = "deg"
  
  
  #unify angles
  data_fp$deg[data_fp$deg==225] = 135
  data_fp$deg[data_fp$deg==315] = 45
  data_fp$deg[data_fp$deg==270] = 90
  #rename variables
  data_fp$reactionTime=data_fp$response_time
  data_fp$response_time=NULL
  
  #####  identify and remove outlier based on rotation reaction time (mean+-3sd)   #####
  #data=sortOutliers(verbose,data,outlierFactor)
  data_fp=sortOutliers(verbose,data_fp,outlierFactor)
  
  
  data_fp$ID = as.factor(data_fp$ID)
  data_fp$block = as.factor(data_fp$block)
  
  data_rpe$block[data_rpe$block=="visualQuestionnaire"] = "visual"
  data_rpe$block[data_rpe$block=="visualStopQuestionnaire"] = "visualStop"
  data_rpe$block[data_rpe$block=="mentalQuestionnaire"] = "mental"
  
  return(list(data, data_fp, data_rpe))
}
#participantNumber 157 (all)
getData=function(data,data_fp,data_rpe,participantNumber){
  data_list = load_and_prep_data(data,data_fp,data_rpe)
  
  data = data_list[[1]]
  data_fp = data_list[[2]]
  data_rpe = data_list[[3]]
  
  # as.numeric(sort(as.numeric(gsub("VP", "", unique(data$ID)))))
  # as.numeric(sort(as.numeric(gsub("VP", "", unique(data_fp$ID)))))
  # 
  # length(as.numeric(sort(as.numeric(gsub("VP", "", unique(data$ID))))))
  # length(as.numeric(sort(as.numeric(gsub("VP", "", unique(data_fp$ID))))))
  # 
  # as.numeric(sort(as.numeric(gsub("VP", "", unique(data_fp$ID))))) %in% as.numeric(sort(as.numeric(gsub("VP", "", unique(data$ID)))))
  # as.numeric(sort(as.numeric(gsub("VP", "", unique(data$ID))))) %in% as.numeric(sort(as.numeric(gsub("VP", "", unique(data_fp$ID)))))
  # 
  # all(as.numeric(sort(as.numeric(gsub("VP", "", unique(data_fp$ID))))) %in% as.numeric(sort(as.numeric(gsub("VP", "", unique(data$ID))))))
  # all(as.numeric(sort(as.numeric(gsub("VP", "", unique(data$ID))))) %in% as.numeric(sort(as.numeric(gsub("VP", "", unique(data_fp$ID))))))
  # 
  #dataset with all participants
  IDs_to_30 = as.character(1:participantNumber)
  data = data[gsub("VP", "", data$ID) %in% IDs_to_30,]
  data_fp = data_fp[gsub("VP", "", data_fp$ID) %in% IDs_to_30,]
  data_rpe = data_rpe[gsub("VP", "", data_rpe$ID) %in% IDs_to_30,]
  
  # length(unique(data$ID))
  # length(unique(data_fp$ID))
  # length(unique(data_rpe$ID))
  return(list(data, data_fp, data_rpe))
}

##### get ids with performance below chance level in data
getBelowChanceIdCond=function(data, verbose=2){
  below_chance_id_cond = c()
  below_chance_id = c()
  below_chance = FALSE
  for (id in unique(data$ID)) {
    for (cond in unique(data$block)) {
      if(verbose>2){
        print(paste(id, cond))
        print(mean(data$correct[data$ID==id & data$block==cond]))
      }
      if(mean(data$correct[data$ID==id & data$block==cond])<0.5){
        if(verbose>2){
          print("performance below chance level!")
        }
        below_chance_id_cond = c(below_chance_id_cond, paste(id, cond))
        below_chance_id = c(below_chance_id, id)
        below_chance = TRUE
      }
      if(verbose>2){
        print("------------------------")
      }
    }
  }
  if(below_chance & verbose >1){
    print("performance below chance level for ID-condition combinations: ")
    print(below_chance_id_cond)
  }else if(verbose>1){
    print("no performance is below chance level.")
  }
  return(below_chance_id)
}

