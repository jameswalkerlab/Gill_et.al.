#BRITTANY LEGER


#PURPOSE----------------------------------------------------------------------------------------------------------------------------------------------------------
#to merge the rhythmicity data in meta_data_analysis by genotype, to produce a percentage of each type subset of rhythmicity score per genotype

#SET UP and DEFINE VARIABLES----------------------------------------------------------------------------------------------------------------------------------------------------------
#setwd
folder <- ""
v<-"23" #version number of analysis code
#"xx" will need to be changed to your proper working directory.
setwd(paste("xx/", folder, "/Data/Analysis (Data)/", subfolder, sep="", collapse=NULL))

#MANUALLY DEFINE FROM META DATA
startSum<-10 #which column do you start doing analysis (aka where do your numerical values start)- set to NFlies column in meta data (discluding x column, index from 1) 


#META_DATA aggregate by genotype for rhythmic analysis MERGE
meta_data <- read.csv("meta_data_analysis.csv",header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
meta_data$X<-NULL
m<-aggregate(meta_data[startSum:ncol(meta_data)], by=list(meta_data$genotype), FUN=sum)

m$Net_Rhythmic_Percent<-(m$Net_Rhythmic/m$NFlies)*100
m$Arrhythmic_Percent<-(m$Arrhythmic/m$NFlies)*100
m$Weak_Rhythmic_Percent<-(m$Weak_Rhythmic/m$NFlies)*100
m$Mid_Rhythmic_Percent<-(m$Mid_Rhythmic/m$NFlies)*100
m$Strong_Rhythmic_Percent<-(m$Strong_Rhythmic/m$NFlies)*100

#change column names to add "genotype"
names(m)[1]<-"genotype"

#WRITE ALL FILES----------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(m, file="rhythmic_analysis_summed.csv")
