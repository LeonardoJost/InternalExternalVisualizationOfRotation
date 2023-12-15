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
}

########## Main Hypothesis ##########
##### Mean Amplitude #####
# fitting various bf-models for comparison
bf_dbintid_dbid = c()
bf_did_dbintid = c()
bf_did_dbid = c()
bf_did_id = c()
bf_did_bid = c()
scaling_factor_fix = 0.5
scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.01, to = 0.11, by = 0.01)

start_time <- Sys.time()
for (scaling_factor_fix in scaling_factors) {
  loop_start_time <- Sys.time()
  
  if(T){
    set.seed(1234)
    bf_full = lmBF(MeanAmplitude ~ deg + block + deg:block + ID,
                   data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_full_no_int = lmBF(MeanAmplitude ~ deg + block + ID,
                          data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_deg = lmBF(MeanAmplitude ~ deg + ID,
                  data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_block = lmBF(MeanAmplitude ~ block + ID,
                    data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_ID = lmBF(MeanAmplitude ~ ID,
                 data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_lm = c(bf_full, bf_full_no_int, bf_deg, bf_block, bf_ID)
    
    
    bf_comp_mat = matrix(NA, nrow = 4, ncol = 4)
    colnames(bf_comp_mat) = row.names(bf_comp_mat) = c("d*b+id", "d+b+id", "d+id", "b+id")
    
    for (i in 1:4) {
      if(i==1){
        #print("log10(BF) of respective model against random intercept-only model:")
      }
      buffer = paste(rep("-", (str_length(unname(names(bf_lm)$numerator[1]))-str_length(unname(names(bf_lm)$numerator[i])))),sep="",collapse="")
      #print(paste0("log10(BF): ", format(round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3),nsmall=3), "   ", ">>", 
       #            unname(names(bf_lm)$numerator[i]), "<<", buffer, " against ", ">>", unname(names(bf_lm)$numerator[5]), "<<"))
      bf_comp_mat[i,i] = round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3)
    }
    for (i in 1:3) {
      for (j in (i+1):4) {
        bf_comp_mat[i,j] = bf_comp_mat[i,i] - bf_comp_mat[j,j]
        bf_comp_mat[j,i] = bf_comp_mat[j,j] - bf_comp_mat[i,i]
      }
    }
    bf_dbintid_dbid = c(bf_dbintid_dbid,bf_comp_mat[2,1])
    bf_did_dbintid = c(bf_did_dbintid,bf_comp_mat[3,1])
    bf_did_dbid = c(bf_did_dbid,bf_comp_mat[3,2])
    bf_did_id = c(bf_did_id,bf_comp_mat[3,3])
    bf_did_bid = c(bf_did_bid,bf_comp_mat[3,4])
    #cat("DV: Mean amplitude \nlog10(BF) against ID-only model in the diagonal. \nModel comparison BFs inupper/lower triangle: row - column")
    #bf_comp_mat
  }
  print(paste("Mean Amplitude analysis for scaling factor", scaling_factor_fix, "done."))
  print(paste("loop duration:", Sys.time()-loop_start_time))
}
print(paste("total duration:", Sys.time()-start_time))


# plot results of robustness region analysis for Mean Amplitude
min_bf = min(c(bf_did_dbintid,bf_did_dbid,bf_did_id,bf_did_bid))
max_bf = max(c(bf_did_dbintid,bf_did_dbid,bf_did_id,bf_did_bid))
plot(scaling_factors, bf_did_dbintid, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
lines(scaling_factors, bf_did_dbid, col='red')
lines(scaling_factors, bf_did_id, col='green')
lines(scaling_factors, bf_did_bid, col='violet')
abline(h=0.477,lty=2)
abline(h=-0.477, lty=2)
abline(v=0.5,lty=2)

#robustness region for interaction
plot(scaling_factors, bf_dbintid_dbid, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")

# Strongest model: d+id log10(BFs) >= 3.301 . Robustness region for d+id > all others: r[0.02,1<]
# Interaction: d*b+id vs d+b+id log10(BF) = 2.558. Robustness region: r[<0.01,1<]

##### Sway Velocity #####
# fitting various bf-models for comparison
bf_dbintid_dbid = c()
bf_id_dbintid = c()
bf_id_dbid = c()
bf_id_did = c()
bf_id_bid = c()

scaling_factor_fix = 0.5
scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.11, to = 0.21, by = 0.01)

start_time <- Sys.time()
for (scaling_factor_fix in scaling_factors) {
  loop_start_time <- Sys.time()
  if(T){
    set.seed(1234)
    bf_full = lmBF(SwayVelocity ~ deg + block + deg:block + ID,
                   data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_full_no_int = lmBF(SwayVelocity ~ deg + block + ID,
                          data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_deg = lmBF(SwayVelocity ~ deg + ID,
                  data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_block = lmBF(SwayVelocity ~ block + ID,
                    data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_ID = lmBF(SwayVelocity ~ ID,
                 data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_lm = c(bf_full, bf_full_no_int, bf_deg, bf_block, bf_ID)
    
    bf_comp_mat = matrix(NA, nrow = 4, ncol = 4)
    colnames(bf_comp_mat) = row.names(bf_comp_mat) = c("d*b+id", "d+b+id", "d+id", "b+id")
    
    for (i in 1:4) {
      if(i==1){
        #print("log10(BF) of respective model against random intercept-only model:")
      }
      buffer = paste(rep("-", (str_length(unname(names(bf_lm)$numerator[1]))-str_length(unname(names(bf_lm)$numerator[i])))),sep="",collapse="")
      #print(paste0("log10(BF): ", format(round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3),nsmall=3), "   ", ">>", 
        #           unname(names(bf_lm)$numerator[i]), "<<", buffer, " against ", ">>", unname(names(bf_lm)$numerator[5]), "<<"))
      bf_comp_mat[i,i] = round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3)
    }
    for (i in 1:3) {
      for (j in (i+1):4) {
        bf_comp_mat[i,j] = bf_comp_mat[i,i] - bf_comp_mat[j,j]
        bf_comp_mat[j,i] = bf_comp_mat[j,j] - bf_comp_mat[i,i]
      }
    }
    bf_dbintid_dbid = c(bf_dbintid_dbid,bf_comp_mat[2,1])
    bf_id_dbintid = c(bf_id_dbintid,bf_comp_mat[1,1])
    bf_id_dbid = c(bf_id_dbid,bf_comp_mat[2,2])
    bf_id_did = c(bf_id_did,bf_comp_mat[3,3])
    bf_id_bid = c(bf_id_bid,bf_comp_mat[4,4])
    #cat("DV: Sway velocity \nlog10(BF) against ID-only model in the diagonal. \nModel comparison BFs inupper/lower triangle: row - column")
    #bf_comp_mat
  }
  print(paste("Sway Velocity analysis for scaling factor", scaling_factor_fix, "done."))
  print(paste("loop duration:", Sys.time()-loop_start_time))
}
print(paste("total duration:", Sys.time()-start_time))

# plot results of robustness region analysis for Sway Velocity
min_bf = min(c(bf_id_dbintid,bf_id_dbid,bf_id_did,bf_id_bid))
max_bf = max(c(bf_id_dbintid,bf_id_dbid,bf_id_did,bf_id_bid))
plot(scaling_factors, bf_id_dbintid, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
lines(scaling_factors, bf_id_dbid, col='red')
lines(scaling_factors, bf_id_did, col='green')
lines(scaling_factors, bf_id_bid, col='violet')
abline(h=0.477,lty=2)
abline(h=-0.477, lty=2)
abline(v=0.5,lty=2)

# Strongest model: id log10(BFs) <= -1.410. Robustness region for id > all others: r[0.17,1<]
# Interaction: d*b+id vs d+b+id log10(BF) = 3.170. Robustness region: r[<0.01,1<]

#### Max Range X ####
# fitting various bf-models for comparison
bf_dbintid_dbid = c()
bf_did_dbintid = c()
bf_did_dbid = c()
bf_did_id = c()
bf_did_bid = c()
scaling_factor_fix = 0.5
scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.01, to = 0.11, by = 0.01)

start_time <- Sys.time()
for (scaling_factor_fix in scaling_factors) {
  loop_start_time <- Sys.time()
  
  if(T){
    set.seed(1234)
    bf_full = lmBF(MaxRangeX ~ deg + block + deg:block + ID,
                   data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_full_no_int = lmBF(MaxRangeX ~ deg + block + ID,
                          data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_deg = lmBF(MaxRangeX ~ deg + ID,
                  data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_block = lmBF(MaxRangeX ~ block + ID,
                    data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_ID = lmBF(MaxRangeX ~ ID,
                 data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_lm = c(bf_full, bf_full_no_int, bf_deg, bf_block, bf_ID)
    
    
    bf_comp_mat = matrix(NA, nrow = 4, ncol = 4)
    colnames(bf_comp_mat) = row.names(bf_comp_mat) = c("d*b+id", "d+b+id", "d+id", "b+id")
    
    for (i in 1:4) {
      if(i==1){
        #print("log10(BF) of respective model against random intercept-only model:")
      }
      buffer = paste(rep("-", (str_length(unname(names(bf_lm)$numerator[1]))-str_length(unname(names(bf_lm)$numerator[i])))),sep="",collapse="")
      #print(paste0("log10(BF): ", format(round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3),nsmall=3), "   ", ">>", 
      #            unname(names(bf_lm)$numerator[i]), "<<", buffer, " against ", ">>", unname(names(bf_lm)$numerator[5]), "<<"))
      bf_comp_mat[i,i] = round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3)
    }
    for (i in 1:3) {
      for (j in (i+1):4) {
        bf_comp_mat[i,j] = bf_comp_mat[i,i] - bf_comp_mat[j,j]
        bf_comp_mat[j,i] = bf_comp_mat[j,j] - bf_comp_mat[i,i]
      }
    }
    bf_dbintid_dbid = c(bf_dbintid_dbid,bf_comp_mat[2,1])
    bf_did_dbintid = c(bf_did_dbintid,bf_comp_mat[3,1])
    bf_did_dbid = c(bf_did_dbid,bf_comp_mat[3,2])
    bf_did_id = c(bf_did_id,bf_comp_mat[3,3])
    bf_did_bid = c(bf_did_bid,bf_comp_mat[3,4])
    #cat("DV: Max range X \nlog10(BF) against ID-only model in the diagonal. \nModel comparison BFs inupper/lower triangle: row - column")
    #bf_comp_mat
  }
  print(paste("Max Range X analysis for scaling factor", scaling_factor_fix, "done."))
  print(paste("loop duration:", Sys.time()-loop_start_time))
}
print(paste("total duration:", Sys.time()-start_time))

# plot results of robustness region analysis for Max Range X
min_bf = min(c(bf_did_dbintid,bf_did_dbid,bf_did_id,bf_did_bid))
max_bf = max(c(bf_did_dbintid,bf_did_dbid,bf_did_id,bf_did_bid))
plot(scaling_factors, bf_did_dbintid, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
lines(scaling_factors, bf_did_dbid, col='red')
lines(scaling_factors, bf_did_id, col='green')
lines(scaling_factors, bf_did_bid, col='violet')
abline(h=0.477,lty=2)
abline(h=-0.477, lty=2)
abline(v=0.5,lty=2)

# Strongest model: d+id log10(BFs) >= 2.279. Robustness region for d+id > all others: r[0.05,1<]
# Interaction: d*b+id vs d+b+id log10(BF) = 2.161. Robustness region: r[<0.01,1<]

#### Max Range Y ####
# fitting various bf-models for comparison
bf_dbintid_dbid = c()
bf_did_dbintid = c()
bf_did_dbid = c()
bf_did_id = c()
bf_did_bid = c()
scaling_factor_fix = 0.5
scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.01, to = 0.11, by = 0.01)

start_time <- Sys.time()
for (scaling_factor_fix in scaling_factors) {
  loop_start_time <- Sys.time()
  
  if(T){
    set.seed(1234)
    bf_full = lmBF(MaxRangeY ~ deg + block + deg:block + ID,
                   data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_full_no_int = lmBF(MaxRangeY ~ deg + block + ID,
                          data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_deg = lmBF(MaxRangeY ~ deg + ID,
                  data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_block = lmBF(MaxRangeY ~ block + ID,
                    data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_ID = lmBF(MaxRangeY ~ ID,
                 data = data_fp, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_lm = c(bf_full, bf_full_no_int, bf_deg, bf_block, bf_ID)
    
    
    bf_comp_mat = matrix(NA, nrow = 4, ncol = 4)
    colnames(bf_comp_mat) = row.names(bf_comp_mat) = c("d*b+id", "d+b+id", "d+id", "b+id")
    
    for (i in 1:4) {
      if(i==1){
        #print("log10(BF) of respective model against random intercept-only model:")
      }
      buffer = paste(rep("-", (str_length(unname(names(bf_lm)$numerator[1]))-str_length(unname(names(bf_lm)$numerator[i])))),sep="",collapse="")
      #print(paste0("log10(BF): ", format(round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3),nsmall=3), "   ", ">>", 
      #            unname(names(bf_lm)$numerator[i]), "<<", buffer, " against ", ">>", unname(names(bf_lm)$numerator[5]), "<<"))
      bf_comp_mat[i,i] = round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3)
    }
    for (i in 1:3) {
      for (j in (i+1):4) {
        bf_comp_mat[i,j] = bf_comp_mat[i,i] - bf_comp_mat[j,j]
        bf_comp_mat[j,i] = bf_comp_mat[j,j] - bf_comp_mat[i,i]
      }
    }
    bf_dbintid_dbid = c(bf_dbintid_dbid,bf_comp_mat[2,1])
    bf_did_dbintid = c(bf_did_dbintid,bf_comp_mat[3,1])
    bf_did_dbid = c(bf_did_dbid,bf_comp_mat[3,2])
    bf_did_id = c(bf_did_id,bf_comp_mat[3,3])
    bf_did_bid = c(bf_did_bid,bf_comp_mat[3,4])
    #cat("DV: Max range Y \nlog10(BF) against ID-only model in the diagonal. \nModel comparison BFs inupper/lower triangle: row - column")
    #bf_comp_mat
  }
  print(paste("Max Range Y analysis for scaling factor", scaling_factor_fix, "done."))
  print(paste("loop duration:", Sys.time()-loop_start_time))
}
print(paste("total duration:", Sys.time()-start_time))

# plot results of robustness region analysis for Max Range Y
min_bf = min(c(bf_did_dbintid,bf_did_dbid,bf_did_id,bf_did_bid))
max_bf = max(c(bf_did_dbintid,bf_did_dbid,bf_did_id,bf_did_bid))
plot(scaling_factors, bf_did_dbintid, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
lines(scaling_factors, bf_did_dbid, col='red')
lines(scaling_factors, bf_did_id, col='green')
lines(scaling_factors, bf_did_bid, col='violet')
abline(h=0.477,lty=2)
abline(h=-0.477, lty=2)
abline(v=0.5,lty=2)

# Strongest model: d+id log10(BFs) >= 3.240 . Robustness region for d+id > all others: r[0.02,1<]
# Interaction: d*b+id vs d+b+id log10(BF) = 2.281. Robustness region: r[<0.01,1<]

########## Secondary Hypotheses ##########
##### Hypothesis S1: Cognitive Effort #####
data_cRPE = data_rpe[data_rpe$type=="cRPE",]

scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.01, to = 0.11, by = 0.01)
scaling_factor_fix = 0.5
bfs_bid_id = c()
for (scaling_factor_fix in scaling_factors) {
  set.seed(1234)
  bf_block = lmBF(value ~ block + ID, data = data_cRPE, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf_ID = lmBF(value ~ ID, data = data_cRPE, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf = bf_block/bf_ID
  bfs_bid_id = c(bfs_bid_id, exp(bf@bayesFactor$bf))
}
min_bf = min(bfs_bid_id)
max_bf = max(bfs_bid_id)
plot(scaling_factors, bfs_bid_id, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
abline(h=1/3,lty=2)
abline(h=3, lty=2)
abline(v=0.5,lty=2)

# b+id model preferred BF = 1.57019e+86. Robustness region for b+id > id: r[<0.01,>1]

# pairwise ANOVA-comparison with random effect ID
# mental vs. visual
ce_ment_vis = data_cRPE[data_cRPE$block!="visualStop",]

scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.01, to = 0.11, by = 0.01)
scaling_factor_fix = 0.5
bfs_bid_id = c()
for (scaling_factor_fix in scaling_factors) {
  set.seed(1234)
  bf_block = lmBF(value ~ block + ID, data = ce_ment_vis, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf_ID = lmBF(value ~ ID, data = ce_ment_vis, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf = bf_block/bf_ID
  bfs_bid_id = c(bfs_bid_id, exp(bf@bayesFactor$bf))
}
min_bf = min(bfs_bid_id)
max_bf = max(bfs_bid_id)
plot(scaling_factors, bfs_bid_id, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
abline(h=1/3,lty=2)
abline(h=3, lty=2)
abline(v=0.5,lty=2)

# b+id model preferred BF = 4.941869e+63. Robustness region for d+id > id: r[<0.01,>1]


# mental vs. visualStop
ce_ment_visStop = data_cRPE[data_cRPE$block!="visual",]

scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.01, to = 0.11, by = 0.01)
scaling_factor_fix = 0.5
bfs_bid_id = c()
for (scaling_factor_fix in scaling_factors) {
  set.seed(1234)
  bf_block = lmBF(value ~ block + ID, data = ce_ment_visStop, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf_ID = lmBF(value ~ ID, data = ce_ment_visStop, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf = bf_block/bf_ID
  bfs_bid_id = c(bfs_bid_id, exp(bf@bayesFactor$bf))
}
min_bf = min(bfs_bid_id)
max_bf = max(bfs_bid_id)
plot(scaling_factors, bfs_bid_id, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
abline(h=1/3,lty=2)
abline(h=3, lty=2)
abline(v=0.5,lty=2)

# b+id model preferred BF = 2.339962e+60. Robustness region for d+id > id: r[<0.01,>1]

# visual vs. visualStop
ce_vis_visStop = data_cRPE[data_cRPE$block!="mental",]

scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.01, to = 0.13, by = 0.01)
scaling_factor_fix = 0.5
bfs_bid_id = c()
for (scaling_factor_fix in scaling_factors) {
  set.seed(1234)
  bf_block = lmBF(value ~ block + ID, data = ce_vis_visStop, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf_ID = lmBF(value ~ ID, data = ce_vis_visStop, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf = bf_block/bf_ID
  bfs_bid_id = c(bfs_bid_id, exp(bf@bayesFactor$bf))
}
min_bf = min(bfs_bid_id)
max_bf = max(bfs_bid_id)
plot(scaling_factors, bfs_bid_id, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
abline(h=1/3,lty=2)
abline(h=3, lty=2)
abline(v=0.5,lty=2)

# id model preferred BF = 0.08961437. Robustness region for id > b+id: r[0.12,>1]

##### Hypothesis S2: Physical Effort #####
data_RPE = data_rpe[data_rpe$type=="RPE",]

scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.51, to = 0.71, by = 0.01)
scaling_factor_fix = 0.5
bfs_bid_id = c()
for (scaling_factor_fix in scaling_factors) {
  set.seed(1234)
  bf_block = lmBF(value ~ block + ID, data = data_RPE, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf_ID = lmBF(value ~ ID, data = data_RPE, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf = bf_block/bf_ID
  bfs_bid_id = c(bfs_bid_id, exp(bf@bayesFactor$bf))
}
min_bf = min(bfs_bid_id)
max_bf = max(bfs_bid_id)
plot(scaling_factors, bfs_bid_id, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
abline(h=1/3,lty=2)
abline(h=3, lty=2)
abline(v=0.5,lty=2)

# Evidence inconclusive BF = 0.5124945. Robustness region for this: r[0.15,0.63]


##### Hypothesis S3: Reaction Time #####
data_rt = data[data$correct==1 & data$outlier==FALSE,]

ggplot(data_rt, aes(x=deg, y=reactionTime, color = block, group = block)) +
  geom_point(stat='summary', fun=mean) +
  stat_summary(fun=mean, geom="line") +
  scale_x_continuous("degree", labels = as.character(c(45,90,135,180)), 
                     breaks = c(45,90,135,180))+
  coord_cartesian(ylim=c(0, 3000))

current_rt_data = data_rt
current_rt_data = data_rt[data_rt$block!="visualStop",]
current_rt_data = data_rt[data_rt$block!="visual",]
current_rt_data = data_rt[data_rt$block!="mental",]

bf_dbintid_did = c()
bf_dbintid_dbid = c()
bf_dbintid_id = c()
bf_dbintid_bid = c()

bf_dbid_did = c()
bf_dbid_bid = c()

scaling_factor_fix = 0.5
scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.01, to = 0.11, by = 0.01)

start_time <- Sys.time()
for (scaling_factor_fix in scaling_factors) {
  loop_start_time <- Sys.time()
  if(T){
    set.seed(1234)
    bf_full = lmBF(reactionTime ~ deg + block + deg:block + ID,
                   data = current_rt_data, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_full_no_int = lmBF(reactionTime ~ deg + block + ID,
                          data = current_rt_data, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_deg = lmBF(reactionTime ~ deg + ID,
                  data = current_rt_data, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_block = lmBF(reactionTime ~ block + ID,
                    data = current_rt_data, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_ID = lmBF(reactionTime ~ ID,
                 data = current_rt_data, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_lm = c(bf_full, bf_full_no_int, bf_deg, bf_block, bf_ID)
    
    bf_comp_mat = matrix(NA, nrow = 4, ncol = 4)
    colnames(bf_comp_mat) = row.names(bf_comp_mat) = c("d*b+id", "d+b+id", "d+id", "b+id")
    
    # Bayes factors are to big (>x*10**308) so R handles them as inf. Calculate relative BFs from output.
    # extract factors an exponents for base 10 from output
    factors = c()
    exponents = c()
    for(i in 1:length(bf_lm)){
      tmp = vector('character')
      con    <- textConnection('tmp', 'wr', local = TRUE)
      sink(con)
      print(bf_lm[i])
      sink()
      close(con)
      extrctd_bf = tail(strsplit(tmp[3], "\\ ")[[1]][-2],2)[1]
      factors = c(factors, strsplit(extrctd_bf, "e\\+")[[1]][1])
      exponents = c(exponents, strsplit(extrctd_bf, "e\\+")[[1]][2])
    }
    # calculate BFs for each model compared to id-only model
    factors = as.numeric(factors)
    exponents = as.numeric(exponents)
    bf_to_id_only_fac = factors[1:4] / factors[5]
    bf_to_id_only_exp = exponents[1:4] - exponents[5]
    
    
    for (i in 1:4) {
      bf_comp_mat[i,i] = bf_to_id_only_fac[i] * 10**bf_to_id_only_exp[i]
      if(i<4){
        for (j in (i+1):4) {
          bf_comp_mat[i,j] = (bf_to_id_only_fac[i]/bf_to_id_only_fac[j]) * 10**(bf_to_id_only_exp[i]-bf_to_id_only_exp[j])
          bf_comp_mat[j,i] = (bf_to_id_only_fac[j]/bf_to_id_only_fac[i]) * 10**(bf_to_id_only_exp[j]-bf_to_id_only_exp[i])
        }
      }
    }
    bf_comp_mat = log10(bf_comp_mat)
    #cat("DV: Reaction time \nlog10(BF) against ID-only model in the diagonal. \nModel comparison BFs inupper/lower triangle: row - column")
    #bf_comp_mat
  }
  bf_dbintid_id = c(bf_dbintid_id,bf_comp_mat[1,1])
  bf_dbintid_dbid = c(bf_dbintid_dbid,bf_comp_mat[1,2])
  bf_dbintid_did = c(bf_dbintid_did,bf_comp_mat[1,3])
  bf_dbintid_bid = c(bf_dbintid_bid,bf_comp_mat[1,4])
  bf_dbid_did = c(bf_dbid_did,bf_comp_mat[2,3])
  bf_dbid_bid = c(bf_dbid_bid,bf_comp_mat[2,4])
  
  
  
  print(paste("Reaction time analysis for scaling factor", scaling_factor_fix, "done."))
  print(paste("loop duration:", Sys.time()-loop_start_time))
}
print(paste("total duration:", Sys.time()-start_time))

plot(scaling_factors, bf_dbintid_dbid, type='l', col='blue', ylab = "BF", xlab = "r")

# Strongest model: d*b+id log10(BFs) >= 2.219. Robustness region for d*b+id > all others: r[<0.01,>1]
# Main Effect condition: d+b+id vs d+id log10(BF) = Inf. Robustness region: r[<0.01,>1]
# Main Effect Degree: d+b+id vs b+id log10(BF) = Inf. Robustness region: r[<0.01,>1]

# Strongest model for comparison mental-visual: d*b+id log10(BFs) >= 2.936. Robustness region for d*b+id > all others: r[<0.01,>1]
# Main Effect condition: d+b+id vs d+id log10(BF) = Inf. Robustness region: r[<0.01,>1]

# Strongest model for comparison mental-visualStop: d*b+id log10(BFs) >= 1.123. Robustness region for d*b+id > all others: r[<0.01,>1]
# Main Effect condition: d+b+id vs d+id log10(BF) = Inf. Robustness region: r[<0.01,>1]

# Strongest model for comparison visual-visualStop: d+b+id log10(BFs) >= 1.110. Robustness region for d+b+id > all others: r[<0.01,>1]
# Interaction: d*b+id vs d+b+id log10(BF) = -1.572 Robustness region: r[<0.01,>1]


##### Hypothesis S4: Pupil Size #####
pupildata = data[data$correct==1 & data$outlier==FALSE,]

trials = c()
for (id in unique(pupildata$ID)) {
  print(id)
  print(paste("total trials: ", length(pupildata$PupilSizeDif[pupildata$ID==id])))
  print(paste("missing values for pupil size: ", sum(is.na(pupildata$PupilSizeDif[pupildata$ID==id]))))
  trials=c(trials, length(pupildata$PupilSizeDif[pupildata$ID==id]))
  print("------------------------------")
}
for (id in unique(pupildata$ID)) {
  if(any(!is.na(pupildata$baselinePupilSize[pupildata$ID==id]))){
    print(paste("create baseline pupil size histogram for ID", id))
    hist(pupildata$baselinePupilSize[pupildata$ID==id], main = paste0("ID = ", id))
  }else{
    print(paste("no values for baseline pupil size for ID", id))
  }
}

pupildata = pupildata[!is.na(pupildata$PupilSizeDif),]

bf_dbid_dbintid = c()
bf_dbid_did = c()
bf_dbid_id = c()
bf_dbid_bid = c()
scaling_factor_fix = 0.5
scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.8, to = 1, by = 0.01)

start_time <- Sys.time()
for (scaling_factor_fix in scaling_factors) {
  loop_start_time <- Sys.time()
  
  if(T){
    set.seed(1234)
    bf_full = lmBF(PupilSizeDif ~ deg + block + deg:block + ID,
                   data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_full_no_int = lmBF(PupilSizeDif ~ deg + block + ID,
                          data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_deg = lmBF(PupilSizeDif ~ deg + ID,
                  data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_block = lmBF(PupilSizeDif ~ block + ID,
                    data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_ID = lmBF(PupilSizeDif ~ ID,
                 data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_lm = c(bf_full, bf_full_no_int, bf_deg, bf_block, bf_ID)
    
    
    bf_comp_mat = matrix(NA, nrow = 4, ncol = 4)
    colnames(bf_comp_mat) = row.names(bf_comp_mat) = c("d*b+id", "d+b+id", "d+id", "b+id")
    
    for (i in 1:4) {
      if(i==1){
        #print("log10(BF) of respective model against random intercept-only model:")
      }
      buffer = paste(rep("-", (str_length(unname(names(bf_lm)$numerator[1]))-str_length(unname(names(bf_lm)$numerator[i])))),sep="",collapse="")
      #print(paste0("log10(BF): ", format(round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3),nsmall=3), "   ", ">>", 
      #            unname(names(bf_lm)$numerator[i]), "<<", buffer, " against ", ">>", unname(names(bf_lm)$numerator[5]), "<<"))
      bf_comp_mat[i,i] = round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3)
    }
    for (i in 1:3) {
      for (j in (i+1):4) {
        bf_comp_mat[i,j] = bf_comp_mat[i,i] - bf_comp_mat[j,j]
        bf_comp_mat[j,i] = bf_comp_mat[j,j] - bf_comp_mat[i,i]
      }
    }
    bf_dbid_dbintid = c(bf_dbid_dbintid,bf_comp_mat[2,1])
    bf_dbid_did = c(bf_dbid_did,bf_comp_mat[2,3])
    bf_dbid_id = c(bf_dbid_id,bf_comp_mat[2,2])
    bf_dbid_bid = c(bf_dbid_bid,bf_comp_mat[2,4])
    #cat("DV: Mean amplitude \nlog10(BF) against ID-only model in the diagonal. \nModel comparison BFs inupper/lower triangle: row - column")
    #bf_comp_mat
  }
  print(paste("Mean pupil size analysis for scaling factor", scaling_factor_fix, "done."))
  print(paste("loop duration:", Sys.time()-loop_start_time))
}
print(paste("total duration:", Sys.time()-start_time))


# plot results of robustness region analysis for Mean Amplitude
min_bf = min(c(bf_dbid_dbintid,bf_dbid_did,bf_dbid_id,bf_dbid_bid))
max_bf = max(c(bf_dbid_dbintid,bf_dbid_did,bf_dbid_id,bf_dbid_bid))
plot(scaling_factors, bf_dbid_dbintid, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
lines(scaling_factors, bf_dbid_did, col='red')
lines(scaling_factors, bf_dbid_id, col='green')
lines(scaling_factors, bf_dbid_bid, col='violet')
abline(h=0.477,lty=2)
abline(h=-0.477, lty=2)
abline(v=0.5,lty=2)

# Strongest model: d+b+id log10(BFs) >= 0.941. Robustness region for d+b+id > all others: r[<0.01,0.85]

# pairwise ANOVA-comparison with random effect ID


current_data = pupildata[pupildata$block!="visualStop",]
current_data = pupildata[pupildata$block!="visual",]
current_data = pupildata[pupildata$block!="mental",]

scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.21, to = 0.31, by = 0.01)
scaling_factor_fix = 0.5
bfs_dbid_did = c()
for (scaling_factor_fix in scaling_factors) {
  set.seed(1234)
  bf_block = lmBF(PupilSizeDif ~ deg + block + ID, data = current_data, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf_ID = lmBF(PupilSizeDif ~ deg + ID, data = current_data, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf = bf_block/bf_ID
  bfs_dbid_did = c(bfs_dbid_did, exp(bf@bayesFactor$bf))
}
min_bf = min(bfs_dbid_did)
max_bf = max(bfs_dbid_did)
plot(scaling_factors, bfs_dbid_did, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
abline(h=1/3,lty=2)
abline(h=3, lty=2)
abline(v=0.5,lty=2)

# Strongest model for comparison mental-visual: d+b+id BF=193.1779. Robustness region for d+b+id > d+id: r[<0.01,>1]
# Strongest model for comparison mental-visualStop: inconclusive BF=0.408. Robustness region for inconclusive: r[0.06,0.61]
# Strongest model for comparison visual-visualStop: d+id BF=0.144. Robustness region for d+id > d+b+id: r[<0.01,>1]


##### Exploratory Analyses #####
##### Max Pupil Size #####
pupildata = data[data$correct==1 & data$outlier==FALSE,]

trials = c()
for (id in unique(pupildata$ID)) {
  print(id)
  print(paste("total trials: ", length(pupildata$maxPupilSizeDif[pupildata$ID==id])))
  print(paste("missing values for pupil size: ", sum(is.na(pupildata$maxPupilSizeDif[pupildata$ID==id]))))
  trials=c(trials, length(pupildata$maxPupilSizeDif[pupildata$ID==id]))
  print("------------------------------")
}

pupildata = pupildata[!is.na(pupildata$maxPupilSizeDif),]

bf_dbid_dbintid = c()
bf_dbid_did = c()
bf_dbid_id = c()
bf_dbid_bid = c()
scaling_factor_fix = 0.5
scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)

start_time <- Sys.time()
for (scaling_factor_fix in scaling_factors) {
  loop_start_time <- Sys.time()
  
  if(T){
    set.seed(1234)
    bf_full = lmBF(maxPupilSizeDif ~ deg + block + deg:block + ID,
                   data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_full_no_int = lmBF(maxPupilSizeDif ~ deg + block + ID,
                          data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_deg = lmBF(maxPupilSizeDif ~ deg + ID,
                  data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_block = lmBF(maxPupilSizeDif ~ block + ID,
                    data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_ID = lmBF(maxPupilSizeDif ~ ID,
                 data = pupildata, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
    bf_lm = c(bf_full, bf_full_no_int, bf_deg, bf_block, bf_ID)
    
    
    bf_comp_mat = matrix(NA, nrow = 4, ncol = 4)
    colnames(bf_comp_mat) = row.names(bf_comp_mat) = c("d*b+id", "d+b+id", "d+id", "b+id")
    
    for (i in 1:4) {
      if(i==1){
        #print("log10(BF) of respective model against random intercept-only model:")
      }
      buffer = paste(rep("-", (str_length(unname(names(bf_lm)$numerator[1]))-str_length(unname(names(bf_lm)$numerator[i])))),sep="",collapse="")
      #print(paste0("log10(BF): ", format(round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3),nsmall=3), "   ", ">>", 
      #            unname(names(bf_lm)$numerator[i]), "<<", buffer, " against ", ">>", unname(names(bf_lm)$numerator[5]), "<<"))
      bf_comp_mat[i,i] = round(log10(extractBF(bf_lm[i]/bf_lm[5])$bf),digits=3)
    }
    for (i in 1:3) {
      for (j in (i+1):4) {
        bf_comp_mat[i,j] = bf_comp_mat[i,i] - bf_comp_mat[j,j]
        bf_comp_mat[j,i] = bf_comp_mat[j,j] - bf_comp_mat[i,i]
      }
    }
    bf_dbid_dbintid = c(bf_dbid_dbintid,bf_comp_mat[2,1])
    bf_dbid_did = c(bf_dbid_did,bf_comp_mat[2,3])
    bf_dbid_id = c(bf_dbid_id,bf_comp_mat[2,2])
    bf_dbid_bid = c(bf_dbid_bid,bf_comp_mat[2,4])
    #cat("DV: Mean amplitude \nlog10(BF) against ID-only model in the diagonal. \nModel comparison BFs inupper/lower triangle: row - column")
    #bf_comp_mat
  }
  print(paste("Max pupil size analysis for scaling factor", scaling_factor_fix, "done."))
  print(paste("loop duration:", Sys.time()-loop_start_time))
}
print(paste("total duration:", Sys.time()-start_time))


# plot results of robustness region analysis for Mean Amplitude
min_bf = min(c(bf_dbid_dbintid,bf_dbid_did,bf_dbid_id,bf_dbid_bid))
max_bf = max(c(bf_dbid_dbintid,bf_dbid_did,bf_dbid_id,bf_dbid_bid))
plot(scaling_factors, bf_dbid_dbintid, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
lines(scaling_factors, bf_dbid_did, col='red')
lines(scaling_factors, bf_dbid_id, col='green')
lines(scaling_factors, bf_dbid_bid, col='violet')
abline(h=0.477,lty=2)
abline(h=-0.477, lty=2)
abline(v=0.5,lty=2)

# Strongest model: d+b+id log10(BFs) >= 2.992. Robustness region for d+b+id > all others: r[<0.01,>1]

# pairwise ANOVA-comparison with random effect ID
current_data = pupildata[pupildata$block!="visualStop",]
current_data = pupildata[pupildata$block!="visual",]
current_data = pupildata[pupildata$block!="mental",]

scaling_factors = c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1)
scaling_factors = seq(from=0.21, to = 0.31, by = 0.01)
scaling_factor_fix = 0.5
bfs_dbid_did = c()
for (scaling_factor_fix in scaling_factors) {
  set.seed(1234)
  bf_block = lmBF(maxPupilSizeDif ~ deg + block + ID, data = current_data, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf_ID = lmBF(maxPupilSizeDif ~ deg + ID, data = current_data, whichRandom = "ID", rscaleFixed = scaling_factor_fix)
  bf = bf_block/bf_ID
  round(log10(extractBF(bf_block/bf_ID)$bf),digits=3)
  bfs_dbid_did = c(bfs_dbid_did, exp(bf@bayesFactor$bf))
}
min_bf = min(bfs_dbid_did)
max_bf = max(bfs_dbid_did)
plot(scaling_factors, bfs_dbid_did, type='l', col='blue', ylim = c(min_bf,max_bf), ylab = "BF", xlab = "r")
abline(h=1/3,lty=2)
abline(h=3, lty=2)
abline(v=0.5,lty=2)

# Strongest model for comparison mental-visual: d+b+id log10BF=7.908. Robustness region for d+b+id > d+id: r[<0.01,>1]
# Strongest model for comparison mental-visualStop: d+b+id log10BF=18.121 Robustness region for d+b+id > d+id: r[<0.01,>1]
# Strongest model for comparison visual-visualStop: inconclusive log10BF=0.264. Robustness region for inconclusive: r[0.31,>1]



##### Gaze data #####
#for separation by conditions
current_data = data
current_data = data[data$block=="visualStop",]
current_data = data[data$block=="visual",]
current_data = data[data$block=="mental",]
##all three figures focused
sum(current_data$gazeOrderUnique %in% c("lrt","rlt","ltr","rtl","tlr","trl"))/sum(current_data$gazeOrderUnique!="")
#correct trials
sum(current_data$gazeOrderUnique %in% c("lrt","rlt","ltr","rtl","tlr","trl") & current_data$correct==1)/sum(current_data$gazeOrderUnique!="" & current_data$correct==1)

##left viewed before right
#sum(current_data$gazeOrderUnique %in% c("l","lt","tl","lr","lrt","ltr","tlr"))#overall 22842
#sum(current_data$gazeOrderUnique %in% c("r","rt","rl","tr","rlt","rtl","trl"))#overall 13880
sum(current_data$gazeOrderUnique %in% c("l","lt","tl","lr","lrt","ltr","tlr"))/
  (sum(current_data$gazeOrderUnique %in% c("l","lt","tl","lr","lrt","ltr","tlr"))+sum(current_data$gazeOrderUnique %in% c("r","rt","rl","tr","rlt","rtl","trl")))
#correct trials
sum(current_data$gazeOrderUnique %in% c("l","lt","tl","lr","lrt","ltr","tlr") & current_data$correct==1)/
  (sum(current_data$gazeOrderUnique %in% c("l","lt","tl","lr","lrt","ltr","tlr") & current_data$correct==1)+
     sum(current_data$gazeOrderUnique %in% c("r","rt","rl","tr","rlt","rtl","trl") & current_data$correct==1))

#focus on correct alternative
sumFocusCorrect=sum(current_data$gazeCWSecondHalfProportions>0.5,na.rm=T)
sumFocusWrong=sum(current_data$gazeCWSecondHalfProportions<0.5,na.rm=T)
suMFocusEqual=sum(current_data$gazeCWSecondHalfProportions==0.5,na.rm=T)
sumFocusCorrect/(sumFocusCorrect+sumFocusWrong+suMFocusEqual)
sumFocusWrong/(sumFocusCorrect+sumFocusWrong+suMFocusEqual)
suMFocusEqual/(sumFocusCorrect+sumFocusWrong+suMFocusEqual)
#correct trials
sumFocusCorrect=sum(current_data$gazeCWSecondHalfProportions>0.5 & current_data$correct==1,na.rm=T)
sumFocusWrong=sum(current_data$gazeCWSecondHalfProportions<0.5 & current_data$correct==1,na.rm=T)
suMFocusEqual=sum(current_data$gazeCWSecondHalfProportions==0.5 & current_data$correct==1,na.rm=T)
sumFocusCorrect/(sumFocusCorrect+sumFocusWrong+suMFocusEqual)

#incorrect trials
sumFocusCorrect=sum(current_data$gazeCWSecondHalfProportions>0.5 & current_data$correct==0,na.rm=T)
sumFocusWrong=sum(current_data$gazeCWSecondHalfProportions<0.5 & current_data$correct==0,na.rm=T)
suMFocusEqual=sum(current_data$gazeCWSecondHalfProportions==0.5 & current_data$correct==0,na.rm=T)
sumFocusCorrect/(sumFocusCorrect+sumFocusWrong+suMFocusEqual)



