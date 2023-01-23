#BRITTANY LEGER


#PURPOSE----------------------------------------------------------------------------------------------------------------------------------------------------------
#purpose: to make average and stdev for tables as ordered by genotype
#input: requires that you define what monitor/row combinations you want to graph
#input: need to define what monitor the controls are that you're going to compare to
#input: requires the products of sleep_analysis_draftN
#uses ggplot
#Analysis_comparison included for DD comparison

#OPEN PACKAGES AND DIRECTORY----------------------------------------------------------------------------
library(plyr)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(psych)

#SET UP and DEFINE VARIABLES----------------------------------------------------------------------------------------------------------------------------------------------------------

folder <- ""
subfolder<-""
v<-"23" #version number of analysis code
#"xx" will need to be changed to your proper working directory.
setwd(paste("xx/", folder, "/Data/Analysis (Data)/", subfolder, sep="", collapse=NULL))

#MANUALLY DEFINE FROM META DATA
x<-c(1,2,5) #manually choose based on your meta data: want monitor, row, genotype


#READ IN FILES----------------------------------------------------------------------------
meta_data <- read.csv("meta_data_analysis.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
Analysis <- read.csv("analysis_output.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
counts_binned <- read.csv("average_counts_half_hour_bin.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
sleep_binned <- read.csv("average_sleep_hour_bin.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
norm_counts_binned <- read.csv("average_counts_half_hour_bin_normalized.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
norm_sleep_binned <- read.csv("average_sleep_hour_bin_normalized.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)

meta_data$X<-NULL
Analysis$X<-NULL
counts_binned$X<-NULL
sleep_binned$X<-NULL
norm_counts_binned$X<-NULL
norm_sleep_binned$X<-NULL

#MERGE FILES TO PUT IN GENOTYPE----------------------------------------------------------------------------
a<-merge(meta_data[x], Analysis,  by= c("monitor", "row"))
cb<-merge(meta_data[x], counts_binned, by= c("monitor", "row"))
sb<-merge(meta_data[x], sleep_binned, by= c("monitor", "row"))
cbn<-merge(meta_data[x], norm_counts_binned, by= c("monitor", "row"))
sbn<-merge(meta_data[x], norm_sleep_binned, by= c("monitor", "row"))


#WRITE AVERAGE AND SD FILES FOR SLEEP----------------------------------------------------------------------------
x<-c(3, 6:ncol(sb))
sb_geno<-sb[x]
row.names(sb_geno)<- sb$Filename
colnames(sb_geno)<-c("genotype", 0:23)

sbAv<-aggregate(sb_geno[2:25], by=list(sb_geno$genotype), FUN=mean)
colnames(sbAv)<-c("genotype", 0:23)
sbT<-t(sbAv)
sbAvT<-as.data.frame(sbT, row.names=sbT[2:24,1])
colnames(sbAvT)<-sbAv$genotype
sbAvT<-sbAvT[-c(1),]
sbAvT$ZT<-c(0:23)
write.csv(sbAvT, file="sleep_binned_averaged.csv")

sbSD<-aggregate(sb_geno[2:25], by=list(sb_geno$genotype), FUN=sd)
colnames(sbSD)<-c("genotype", 0:23)
sbT<-t(sbSD)
sbSDT<-as.data.frame(sbT, row.names=sbT[2:24,1])
colnames(sbSDT)<-sbSD$genotype
sbSDT<-sbSDT[-c(1),]
sbSDT$ZT<-c(0:23)
write.csv(sbSDT, file="sleep_binned_stdev.csv")

#WRITE AVERAGE AND SD FILES FOR NORMALIZED SLEEP----------------------------------------------------------------------------
x<-c(3, 6:ncol(sbn))
sbn_geno<-sbn[x]
row.names(sbn_geno)<- sbn$Filename
colnames(sbn_geno)<-c("genotype", 0:23)

sbnAv<-aggregate(sbn_geno[2:25], by=list(sbn_geno$genotype), FUN=mean)
colnames(sbnAv)<-c("genotype", 0:23)
sbnT<-t(sbnAv)
sbnAvT<-as.data.frame(sbnT, row.names=sbnT[2:24,1])
colnames(sbnAvT)<-sbnAv$genotype
sbnAvT<-sbnAvT[-c(1),]
sbnAvT$ZT<-c(0:23)
write.csv(sbnAvT, file="sleep_binned_normalized_averaged.csv")

sbnSD<-aggregate(sbn_geno[2:25], by=list(sbn_geno$genotype), FUN=sd)
colnames(sbnSD)<-c("genotype", 0:23)
sbnT<-t(sbnSD)
sbnSDT<-as.data.frame(sbnT, row.names=sbnT[2:24,1])
colnames(sbnSDT)<-sbnSD$genotype
sbnSDT<-sbnSDT[-c(1),]
sbnSDT$ZT<-c(0:23)
write.csv(sbnSDT, file="sleep_binned_normalized_stdev.csv")
#WRITE AVERAGE AND SD FILES FOR COUNTS----------------------------------------------------------------------------
n<-list("genotype","0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5", "12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5")
x<-c(3, 6:ncol(cb))
cb_geno<-cb[x]
row.names(cb_geno)<- cb$Filename
colnames(cb_geno)<-c(n)

cbAv<-aggregate(cb_geno[2:ncol(cb_geno)], by=list(cb_geno$genotype), FUN=mean)
colnames(cbAv)<-n
cbT<-t(cbAv)
cbAvT<-as.data.frame(cbT, row.names=cbT[2:(ncol(cb_geno)-1),1])
colnames(cbAvT)<-cbAv$genotype
cbAvT<-cbAvT[-c(1),]
cbAvT$ZT<-c("0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5", "12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5")
write.csv(cbAvT, file="counts_binned_averaged.csv")

cbSD<-aggregate(cb_geno[2:ncol(cb_geno)], by=list(cb_geno$genotype), FUN=sd)
colnames(cbSD)<-n
cbT<-t(cbSD)
cbSDT<-as.data.frame(cbT, row.names=cbT[2:(ncol(cb_geno)-1),1])
colnames(cbSDT)<-cbSD$genotype
cbSDT<-cbSDT[-c(1),]
cbSDT$ZT<-c("0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5", "12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5")
write.csv(cbSDT, file="counts_binned_stdev.csv")

#WRITE AVERAGE AND SD FILES FOR NORMALIZED COUNTS----------------------------------------------------------------------------
n<-list("genotype","0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5", "12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5")
x<-c(3, 6:ncol(cbn))
cbn_geno<-cbn[x]
row.names(cbn_geno)<- cbn$Filename
colnames(cbn_geno)<-c(n)

cbnAv<-aggregate(cbn_geno[2:ncol(cbn_geno)], by=list(cbn_geno$genotype), FUN=mean)
colnames(cbnAv)<-n
cbnT<-t(cbnAv)
cbnAvT<-as.data.frame(cbnT, row.names=cbnT[2:(ncol(cbn_geno)-1),1])
colnames(cbnAvT)<-cbnAv$genotype
cbnAvT<-cbnAvT[-c(1),]
cbnAvT$ZT<-c("0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5", "12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5")
write.csv(cbnAvT, file="counts_binned_normalized_averaged.csv")

cbnSD<-aggregate(cbn_geno[2:ncol(cbn_geno)], by=list(cbn_geno$genotype), FUN=sd)
colnames(cbnSD)<-n
cbnT<-t(cbnSD)
cbnSDT<-as.data.frame(cbnT, row.names=cbnT[2:(ncol(cbn_geno)-1),1])
colnames(cbnSDT)<-cbnSD$genotype
cbnSDT<-cbnSDT[-c(1),]
cbnSDT$ZT<-c("0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5", "12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5")
write.csv(cbnSDT, file="counts_binned_normalized_stdev.csv")

#WRITE AVERAGE AND SD FILES FOR ANALYSIS_OUTPUT----------------------------------------------------------------------------
n<-as.list(colnames(a))
x<-c(3, 6:ncol(a))
n<-n[x]
a_geno<-a[x]
row.names(a_geno)<- a$Filename
colnames(a_geno)<-c(n)

aAv<-aggregate(a_geno[2:ncol(a_geno)], by=list(a_geno$genotype), FUN=mean)
colnames(aAv)<-n
aT<-t(aAv)
aAvT<-as.data.frame(aT, row.names=aT[2:(ncol(a_geno)-1),1])
colnames(aAvT)<-aAv$genotype
aAvT<-aAvT[-c(1),]
rownames(aAvT)<-colnames(aAv)[2:ncol(aAv)]
write.csv(aAvT, file="analysis_output_averaged.csv")

aSD<-aggregate(a_geno[2:ncol(a_geno)], by=list(a_geno$genotype), FUN=sd)
colnames(aSD)<-n
aT<-t(aSD)
aSDT<-as.data.frame(aT, row.names=aT[2:(ncol(a_geno)-1),1])
colnames(aSDT)<-aSD$genotype
aSDT<-aSDT[-c(1),]
rownames(aSDT)<-colnames(aSD)[2:ncol(aSD)]
write.csv(aSDT, file="analysis_output_stdev.csv")

