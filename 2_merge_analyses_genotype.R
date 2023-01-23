#BRITTANY LEGER


#PURPOSE----------------------------------------------------------------------------------------------------------------------------------------------------------
#FOR MAKING MERGED INDIVIDUAL FLIES WITH LINE NAME- FOR PLOTTING BY GENOTYPE
#ANALYSIS_OUTPUT, SLEEP, SLEEP NORMALIZED, ACTIVITY COUNTS, ACTIVITY COUNTS NORMALIZED


#SET UP and DEFINE VARIABLES----------------------------------------------------------------------------------------------------------------------------------------------------------
#setwd
folder <- ""
subfolder<-""
v<-"23" #version number of analysis code

#"xx" will need to be changed to your proper working directory.
setwd(paste("xx/", folder, "/Data/Analysis (Data)/", subfolder, sep="", collapse=NULL))

#MANUALLY DEFINE FROM META DATA
x<-c(1,2,5) #manually choose based on your meta data: want monitor, row, genotype
#Set as LD or DD
analysis_type<- "XX"
NDays<-4

meta_data <- read.csv("meta_data_analysis.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)


#WRITE ALL FILES----------------------------------------------------------------------------------------------------------------------------------------------------------------

#ANALYSIS_OUTPUT MERGE
Analysis <- read.csv("analysis_output.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
meta_data$X<-NULL
Analysis$X<-NULL
a<-merge(meta_data[x], Analysis,  by= c("monitor", "row"))
a$TrikFolder<-folder
a$AnalysisType<-analysis_type
a$NDays<-NDays
write.csv(a, file="analysis_output_meta_merged.csv")

#SLEEP TRACE MERGE
Analysis <- read.csv("average_sleep_hour_bin.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
meta_data$X<-NULL
Analysis$X<-NULL
a<-merge(meta_data[x], Analysis,  by= c("monitor", "row"))
a$TrikFolder<-folder
a$AnalysisType<-analysis_type
a$NDays<-NDays
write.csv(a, file="sleep_meta_merged.csv")

#NORMALIZED TRACE SLEEP MERGE
Analysis <- read.csv("average_sleep_hour_bin_normalized.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
meta_data$X<-NULL
Analysis$X<-NULL
a<-merge(meta_data[x], Analysis,  by= c("monitor", "row"))
a$TrikFolder<-folder
a$AnalysisType<-analysis_type
a$NDays<-NDays
write.csv(a, file="sleep_norm_meta_merged.csv")

#ACTIVITY TRACE MERGED
Analysis <- read.csv("average_counts_half_hour_bin.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
meta_data$X<-NULL
Analysis$X<-NULL
a<-merge(meta_data[x], Analysis,  by= c("monitor", "row"))
a$TrikFolder<-folder
a$AnalysisType<-analysis_type
a$NDays<-NDays
write.csv(a, file="counts_meta_merged.csv")

#NORMALIZED ACTIVITY TRACE MERGED
Analysis <- read.csv("average_counts_half_hour_bin_normalized.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
meta_data$X<-NULL
Analysis$X<-NULL
a<-merge(meta_data[x], Analysis,  by= c("monitor", "row"))
a$TrikFolder<-folder
a$AnalysisType<-analysis_type
a$NDays<-NDays
write.csv(a, file="counts_norm_meta_merged.csv")
