### Script for plotting data
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

library(ggplot2)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#read data
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

#add stagger for plotting
plotdata = data
stagger = 3
plotdata$deg_stagger[plotdata$block=="mental"] = plotdata$deg[plotdata$block=="mental"] + stagger
plotdata$deg_stagger[plotdata$block=="visual"] = plotdata$deg[plotdata$block=="visual"] - stagger
plotdata$deg_stagger[plotdata$block=="visualStop"] = plotdata$deg[plotdata$block=="visualStop"]

plotdata_fp = data_fp
plotdata_fp$deg_stagger[plotdata_fp$block=="mental"] = plotdata_fp$deg[plotdata_fp$block=="mental"] + stagger
plotdata_fp$deg_stagger[plotdata_fp$block=="visual"] = plotdata_fp$deg[plotdata_fp$block=="visual"] - stagger
plotdata_fp$deg_stagger[plotdata_fp$block=="visualStop"] = plotdata_fp$deg[plotdata_fp$block=="visualStop"]

######## Pretty plotting ########
### Mean amplitude
short_long_fp = summarySE(plotdata_fp, measurevar="MeanAmplitude", groupvars=c("block","deg"))
names(short_long_fp)[names(short_long_fp)=="block"] = "Condition"
short_long_fp$Condition[short_long_fp$Condition=="mental"] = "Mental rotation"
short_long_fp$Condition[short_long_fp$Condition=="visual"] = "Visual rotation"
short_long_fp$Condition[short_long_fp$Condition=="visualStop"] = "Visual rotation with stop"

pd <- position_dodge(8)

from = 2 # set limits and labels for y-axis
to = 2.8
by = 0.1

tiff("figs\\MA.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long_fp, aes(x=deg, y=MeanAmplitude, colour=Condition, group = Condition, linetype=Condition)) + 
  #geom_errorbar(aes(ymin=MeanAmplitude-se, ymax=MeanAmplitude+se), width=10, position = pd, linetype="solid", size=1) +
  geom_line(size=1.1, position = pd) +
  geom_point(position = pd, size=3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
  scale_y_continuous("Mean amplitude (mm)", labels = as.character(seq(from=from,to=to,by=by)), 
                     breaks = seq(from=from,to=to,by=by))+
  scale_linetype_manual(values=c(1,2,3))+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(from, to))
dev.off()

### Sway velocity
short_long_fp = summarySE(plotdata_fp, measurevar="SwayVelocity", groupvars=c("block","deg"))
names(short_long_fp)[names(short_long_fp)=="block"] = "Condition"
short_long_fp$Condition[short_long_fp$Condition=="mental"] = "Mental rotation"
short_long_fp$Condition[short_long_fp$Condition=="visual"] = "Visual rotation"
short_long_fp$Condition[short_long_fp$Condition=="visualStop"] = "Visual rotation with stop"
from = 400 # set limits and labels for y-axis
to = 550
by = 10
tiff("figs\\SV.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long_fp, aes(x=deg, y=SwayVelocity, colour=Condition, group = Condition, linetype=Condition)) + 
  geom_errorbar(aes(ymin=SwayVelocity-se, ymax=SwayVelocity+se), width=10, position = pd, linetype="solid", size=1) +
  geom_line(size=1.1, position = pd) +
  geom_point(position = pd, size=3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
  scale_y_continuous("Sway velocity (mm/s)", labels = as.character(seq(from=from,to=to,by=by)), 
                     breaks = seq(from=from,to=to,by=by))+
  scale_linetype_manual(values=c(1,2,3))+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(from, to))
dev.off()

### Max range X
short_long_fp = summarySE(plotdata_fp, measurevar="MaxRangeX", groupvars=c("block","deg"))
names(short_long_fp)[names(short_long_fp)=="block"] = "Condition"
short_long_fp$Condition[short_long_fp$Condition=="mental"] = "Mental rotation"
short_long_fp$Condition[short_long_fp$Condition=="visual"] = "Visual rotation"
short_long_fp$Condition[short_long_fp$Condition=="visualStop"] = "Visual rotation with stop"
from = 5 # set limits and labels for y-axis
to = 8
by = 1
tiff("figs\\MRX.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long_fp, aes(x=deg, y=MaxRangeX, colour=Condition, group = Condition, linetype=Condition)) + 
  geom_errorbar(aes(ymin=MaxRangeX-se, ymax=MaxRangeX+se), width=10, position = pd, linetype="solid", size=1) +
  geom_line(size=1.1, position = pd) +
  geom_point(position = pd, size=3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
  scale_y_continuous("Max range X (mm)", labels = as.character(seq(from=from,to=to,by=by)), 
                     breaks = seq(from=from,to=to,by=by))+
  scale_linetype_manual(values=c(1,2,3))+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(from, to))
dev.off()

### Max range Y
short_long_fp = summarySE(plotdata_fp, measurevar="MaxRangeY", groupvars=c("block","deg"))
names(short_long_fp)[names(short_long_fp)=="block"] = "Condition"
short_long_fp$Condition[short_long_fp$Condition=="mental"] = "Mental rotation"
short_long_fp$Condition[short_long_fp$Condition=="visual"] = "Visual rotation"
short_long_fp$Condition[short_long_fp$Condition=="visualStop"] = "Visual rotation with stop"
from = 6 # set limits and labels for y-axis
to = 10
by = 1
tiff("figs\\MRY.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long_fp, aes(x=deg, y=MaxRangeY, colour=Condition, group = Condition, linetype=Condition)) + 
  geom_errorbar(aes(ymin=MaxRangeY-se, ymax=MaxRangeY+se), width=10, position = pd, linetype="solid", size=1) +
  geom_line(size=1.1, position = pd) +
  geom_point(position = pd, size=3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
  scale_y_continuous("Max range Y (mm)", labels = as.character(seq(from=from,to=to,by=by)), 
                     breaks = seq(from=from,to=to,by=by))+
  scale_linetype_manual(values=c(1,2,3))+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(from, to))
dev.off()

### subjective effort
short_long = summarySE(data_rpe, measurevar="value", groupvars=c("type","block"))
names(short_long)[names(short_long)=="block"] = "Condition"
short_long$Condition[short_long$Condition=="mental"] = "Mental rotation"
short_long$Condition[short_long$Condition=="visual"] = "Visual rotation"
short_long$Condition[short_long$Condition=="visualStop"] = "Visual rotation with stop"
short_long$type[short_long$type=="cRPE"] = "cognitive"
short_long$type[short_long$type=="RPE"] = "physical"


pd <- position_dodge(8)
from = 1 # set limits and labels for y-axis
to = 5
by = 1
#setwd(file.path("C:", "Users", "Markus", "Nextcloud", "ArbeitUniversitaet",
#                "Leo", "Paper Rotation Postural Stability", "Poster"))
tiff("figs\\rpe.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long, aes(x=type,y=value, fill = Condition)) + 
  #geom_errorbar(aes(ymin=value-se, ymax=value+se), width=10, position = pd, size=1) +
  geom_bar(stat="identity",position = position_dodge()) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  scale_y_continuous("RPE", labels = as.character(seq(from=from,to=to,by=by)), 
                     breaks = seq(from=from,to=to,by=by))+
  xlab("")+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(from, to))
dev.off()

### Reaction time (fp-data)
# short_long_fp = summarySE(plotdata_fp, measurevar="response_time", groupvars=c("block","deg"))
# names(short_long_fp)[names(short_long_fp)=="block"] = "Condition"
# short_long_fp$Condition[short_long_fp$Condition=="mental"] = "Mental rotation"
# short_long_fp$Condition[short_long_fp$Condition=="visual"] = "Visual rotation"
# short_long_fp$Condition[short_long_fp$Condition=="visualStop"] = "Visual rotation with stop"
# from = 1500 # set limits and labels for y-axis
# to = 3000
# by = 500
# tiff("figs\\response_time.tiff", units="in", width=8, height=5, res=300)
# ggplot(short_long_fp, aes(x=deg, y=response_time, colour=Condition, group = Condition, linetype=Condition)) + 
#   geom_errorbar(aes(ymin=response_time-se, ymax=response_time+se), width=10, position = pd, linetype="solid", size=1) +
#   geom_line(size=1.1, position = pd) +
#   geom_point(position = pd, size=3) +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                      axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
#                      axis.title.x = element_text(margin = margin(t = 10)),
#                      axis.title.y = element_text(margin = margin(r = 10)))+
#   scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
#   scale_y_continuous("Reaction time (ms)", labels = as.character(seq(from=from,to=to,by=by)), 
#                      breaks = seq(from=from,to=to,by=by))+
#   coord_cartesian(ylim=c(from, to))
# dev.off()

### Reaction time from leo's data
short_long = summarySE(plotdata, measurevar="reactionTime", groupvars=c("block","deg"))
names(short_long)[names(short_long)=="block"] = "Condition"
short_long$Condition[short_long$Condition=="mental"] = "Mental rotation"
short_long$Condition[short_long$Condition=="visual"] = "Visual rotation"
short_long$Condition[short_long$Condition=="visualStop"] = "Visual rotation with stop"


pd <- position_dodge(8)
from = 1500 # set limits and labels for y-axis
to = 3000
by = 500
#setwd(file.path("C:", "Users", "Markus", "Nextcloud", "ArbeitUniversitaet",
#                "Leo", "Paper Rotation Postural Stability", "Poster"))
tiff("figs\\rt.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long, aes(x=deg, y=reactionTime, colour=Condition, group = Condition, linetype=Condition)) + 
  geom_errorbar(aes(ymin=reactionTime-se, ymax=reactionTime+se), width=10, position = pd, linetype="solid", size=1) +
  geom_line(size=1.1, position = pd) +
  geom_point(position = pd, size=3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
  scale_y_continuous("Reaction time (ms)", labels = as.character(seq(from=from,to=to,by=by)), 
                     breaks = seq(from=from,to=to,by=by))+
  scale_linetype_manual(values=c(1,2,3))+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(from, to))
dev.off()

# rt-diagram starting from 0 on y-axis to paste together with the scaled one
tiff("figs\\rt0.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long, aes(x=deg, y=reactionTime, colour=Condition, group = Condition, linetype=Condition)) + 
  geom_errorbar(aes(ymin=reactionTime-se, ymax=reactionTime+se), width=10, position = pd, linetype="solid", size=1) +
  geom_line(size=1.1, position = pd) +
  geom_point(position = pd, size=3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  theme(legend.position = "none")+
  scale_linetype_manual(values=c(1,2,3))+
  scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
  scale_y_continuous("Reaction time (ms)", labels = as.character(seq(from=0,to=3000,by=500)), 
                     breaks = seq(from=0,to=3000,by=500), limits=c(0,3100))
dev.off()

### Pupil size
short_long = summarySE(plotdata, measurevar="PupilSizeDif", groupvars=c("block","deg"),na.rm=T)
names(short_long)[names(short_long)=="block"] = "Condition"
short_long$Condition[short_long$Condition=="mental"] = "Mental rotation"
short_long$Condition[short_long$Condition=="visual"] = "Visual rotation"
short_long$Condition[short_long$Condition=="visualStop"] = "Visual rotation with stop"


pd <- position_dodge(8)
from = -0.08 # set limits and labels for y-axis
to = -0.05
by = 0.005

tiff("figs\\pupilSizeDif.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long, aes(x=deg, y=PupilSizeDif, colour=Condition, group = Condition, linetype=Condition)) + 
  geom_errorbar(aes(ymin=PupilSizeDif-se, ymax=PupilSizeDif+se), width=10, position = pd, linetype="solid", size=1) +
  geom_line(size=1.1, position = pd) +
  geom_point(position = pd, size=3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
  scale_y_continuous("Pupil size (mean-baseline)", labels = as.character(seq(from=from,to=to,by=by)), 
                     breaks = seq(from=from,to=to,by=by))+
  scale_linetype_manual(values=c(1,2,3))+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(from, to))
dev.off()

### Max pupil size
short_long = summarySE(plotdata, measurevar="maxPupilSizeDif", groupvars=c("block","deg"),na.rm=T)
names(short_long)[names(short_long)=="block"] = "Condition"
short_long$Condition[short_long$Condition=="mental"] = "Mental rotation"
short_long$Condition[short_long$Condition=="visual"] = "Visual rotation"
short_long$Condition[short_long$Condition=="visualStop"] = "Visual rotation with stop"


pd <- position_dodge(8)
from = 0.05 # set limits and labels for y-axis
to = 0.082
by = 0.01

tiff("figs\\pupilSizeDifMax.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long, aes(x=deg, y=maxPupilSizeDif, colour=Condition, group = Condition, linetype=Condition)) + 
  geom_errorbar(aes(ymin=maxPupilSizeDif-se, ymax=maxPupilSizeDif+se), width=10, position = pd, linetype="solid", size=1) +
  geom_line(size=1.1, position = pd) +
  geom_point(position = pd, size=3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
  scale_y_continuous("Pupil size (max-baseline)", labels = as.character(seq(from=from,to=to,by=by)), 
                     breaks = seq(from=from,to=to,by=by))+
  scale_linetype_manual(values=c(1,2,3))+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(from, to))
dev.off()

#histogram of baseline pupil sizes
tiff("figs\\pupil sizes.tiff", units="in", width=8, height=5, res=300)
ggplot(plotdata,aes(baselinePupilSize))+geom_histogram(binwidth=0.5) +xlab("baseline pupil size") + ylab("Count") + theme_bw()
dev.off()

###gaze data
#gaze order
tiff("figs\\gazeOrder.tiff", units="in", width=8, height=5, res=300)
#order levels
plotdata$gazeOrderUnique=factor(plotdata$gazeOrderUnique,levels=c("l","r","lr","rl","lt","rt","tl","tr","lrt","rlt","ltr","rtl","tlr","trl","t",""))
ggplot(plotdata,aes(gazeOrderUnique,fill=type))+geom_bar() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  theme(legend.position = "none")+
  xlab("gaze order") + ylab("Count")
dev.off()

#correct wrong proportions
tiff("figs\\gazeProportionCW.tiff", units="in", width=8, height=5, res=300)
ggplot(plotdata,aes(gazeCWSecondHalfProportions,fill=type))+geom_histogram(binwidth=0.1) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  theme(legend.position = "none")+
  xlab("proportion of gaze on correct answer") + ylab("Count") 
dev.off()

###accuracy
short_long = summarySE(plotdata, measurevar="correct", groupvars=c("block","deg"),na.rm=T)
names(short_long)[names(short_long)=="block"] = "Condition"
short_long$Condition[short_long$Condition=="mental"] = "Mental rotation"
short_long$Condition[short_long$Condition=="visual"] = "Visual rotation"
short_long$Condition[short_long$Condition=="visualStop"] = "Visual rotation with stop"

pd <- position_dodge(8)
from = 0.8 # set limits and labels for y-axis
to = 0.95
by = 0.05

tiff("figs\\acc.tiff", units="in", width=8, height=5, res=300)
ggplot(short_long, aes(x=deg, y=correct, colour=Condition, group = Condition, linetype=Condition)) + 
  geom_errorbar(aes(ymin=correct-se, ymax=correct+se), width=10, position = pd, linetype="solid", size=1) +
  geom_line(size=1.1, position = pd) +
  geom_point(position = pd, size=3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
                     axis.title.x = element_text(margin = margin(t = 10)),
                     axis.title.y = element_text(margin = margin(r = 10)))+
  scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180)) +
  scale_y_continuous("Proportion of correct items", labels = as.character(seq(from=from,to=to,by=by)), 
                     breaks = seq(from=from,to=to,by=by))+
  scale_linetype_manual(values=c(1,2,3))+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(from, to))

# ggplot(plotdata,aes(y=correct,x=deg, fill=block, shape=block,color=block,linetype=block)) + 
#   stat_summary(na.rm=TRUE, fun=mean, geom="line") +
#   stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
#   stat_summary(fun.data=mean_se,geom="errorbar",aes(linetype=NULL)) +
#   labs(x="Angular Disparity",y="Proportion of correct items",color="Condition", linetype="Condition",shape="Condition",fill="Condition")+
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                      axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0),
#                      axis.title.x = element_text(margin = margin(t = 10)),
#                      axis.title.y = element_text(margin = margin(r = 10)))+
#   scale_x_continuous("Rotation angle", labels = c("45","90","135","180"), breaks = c(45,90,135,180))
dev.off()


#combine multiple images into one
combineImages=function(imagesList,rows,columns,outputFile){
  library(magick)
  initImage=image_read(imagesList[1])
  #each row contains columns images
  imageRows=rep(initImage,columns)
  imageColumns=rep(initImage,rows)
  #combine images in rows and columns
  counter=1
  for(i in 1:rows){
    for(j in 1:columns){
      imageRows[j]=image_read(imagesList[counter])
      counter=counter+1
    }
    #append horizontally
    imageRow=image_append(imageRows)
    imageColumns[i]=imageRow
  }
  #append vertically
  image=image_append(imageColumns,stack=TRUE)
  #save
  image_write(image, path = outputFile)
  gc()
}
#add legend centered at bottom
addLegend=function(imageFile, legendFile,outputFile){
  library(magick)
  #read images
  image=image_read(imageFile)
  legend=image_read(legendFile)
  
  legend=image_trim(legend)
  widthLegend=image_info(legend)$width
  imageWidth=image_info(image)$width
  dif=(imageWidth-widthLegend)/2
  legend=image_border(legend,"white",dif)
  image=image_append(c(image,legend),stack=T)
  #save
  image_write(image, path = outputFile)
  gc()
}
#combine specific images
#behavioral measures
combineImages(c("figs/acc.tiff","figs/rt.tiff"),1,2,"figs/behavioral.tiff")
addLegend("figs/behavioral.tiff","figs/conditionlegend.tiff","figs/behvioralLegend.tiff")
#fp measures
combineImages(c("figs/MA.tiff","figs/SV.tiff","figs/MRX.tiff","figs/MRY.tiff"),2,2,"figs/fp.tiff")
addLegend("figs/fp.tiff","figs/conditionlegend.tiff","figs/fpLegend.tiff")
#effort measures
combineImages(c("figs/rpe.tiff","figs/pupilSizeDif.tiff","figs/pupilSizeDifMax.tiff"),1,3,"figs/effort.tiff")
addLegend("figs/effort.tiff","figs/conditionlegend.tiff","figs/effortLegend.tiff")
#gaze patterns
combineImages(c("figs/gazeOrder.tiff","figs/gazeProportionCW.tiff"),1,2,"figs/gaze.tiff")
addLegend("figs/gaze.tiff","figs/typelegend.tiff","figs/gazeLegend.tiff")
