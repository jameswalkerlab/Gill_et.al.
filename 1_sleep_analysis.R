#BRITTANY LEGER

#PURPOSE----------------------------------------------------------------------------------------------------------------------------------------------------------
#the PRIMARY sleep data analysis and consolidation code. MUST BE RUN FIRST.

#OUTPUTS----------------------------------------------------------------------------------------------------------------------------------------------------------
#"date and version.csv"   table of who analyzed it with which version of code on what day
#"analysis_output.csv"    table with values for all sleep measures for individual flies by the location in monitor
#"meta_data_analysis.csv"   table with meta data, as well as rhythmicity raw numbers
#"average_sleep_hour_bin.csv"   table of average min sleep per 1 hour averaged over N days
#"average_sleep_hour_bin_normalized.csv"    table of average min sleep/total min sleep per 1 hour averaged over N days
#"average_counts_half_hour_bin.csv"   table of average number of activity counts per 0.5 hour averaged over N days
#"average_counts_half_hour_bin_normalized.csv"  table of average number of activity counts/total activity counts per 0.5 hour averaged over N days
#all in Zeitgeber time dependent on as defined "lightsON" variable

#NOTES----------------------------------------------------------------------------------------------------------------------------------------------------------
#this code glitches if your monitors are too non-contiguous- I have only had this problem once and I don't know how to fix it.
#this analysis includes everything up through: FFT	PerLength	per_rhythmic	RhythmScore	totalMinSleepPerDay	daySleepAv	nightSleepAv	totalCounts	HammingDistance	NBinPerDay	FragWake	dayCounts	dayBins	dayFrag	nightCounts	nightBins	nightBinLength	nightFrag	duskCounts	dawnCounts	UrnalIndex_Counts	circadianIndex_Counts	BinLength_Min	dayBinLength percent_rhythmic
#including any circadian analysis (period shift and cosinor analysis) requested by Hanna for the Usf project

#SET UP----------------------------------------------------------------------------------------------------------------------------------------------------------
#close all libraries ect. and open all necessary libraries
closeAllConnections()

library(stringi)
library(tidyverse)
library(dplyr)
library(e1071) 


#VARIABLES THAT NEED TO BE DEFINED IN EACH ANALYSIS----------------------------------------------------------------------------------------------------------------------------------------------------------
#set working directory to the location of your data
#set working directory to site of specific data that's being analyzed. All data should be .csv format, and have same naming convention, so all code will work for data that's been properly exported from clocklab. getwd for checking if correct directory
#only folder and subfolder should be manually changed each time it's run
#aka it should be "~/Dropbox (Partners HealthCare)/Trikinetics - SG JW Shared/FOLDER NAME TO BE CHANGED/Data/Analysis (Data)/subfolder"

#example ("/Users/work/Dropbox (Partners HealthCare)/Trikinetics - SG JW Shared/20190913/incubator 1/Data/Analysis (Data)/LD 5 day)
folder <- ""
#folder= trikinetics folder in which Data/Analysis (Data) exists
subfolder<-""
#if you are analyzing this data in multiple different ways, then create a subfolder for each type that is clearly named. This is said folder. If not, leave "" completely empty.
#"xx" will need to be changed to your proper working directory.
setwd(paste("xx/", folder, "/Data/Analysis (Data)/", subfolder, sep="", collapse=NULL))

#set the analysis date, to write into "date and version" file in R_analysis_N folder
dateOfAnalysis<-as.character(Sys.Date())
v<-"23" #version number of analysis code
analysisCodeVersion<-paste("sleep_analysis_draft",v,sep="")
analyst<-"XX" #change to whomever is running this analysis
analysis_meta<-c(dateOfAnalysis, analysisCodeVersion, analyst) 

#folder where the analysis data will be put (within the working directory)
#running an new analysis with the same script as the old analysis will write over the old files (unless you change the folder name for the old analysis manually)
analysis_folder<-paste("R_Analysis",v,sep="_")

#must be defined: for how many days data is being analyzed for
NDays <- 5

#set the hour of lights on and lights off in count table data: military time- you can determine the correct time based on actograms if you're not sure for your data
#manually defined because sometimes lights data is not collected if the wire for the monitor is breaking
lightsON <- 7
lightsOFF <- 12+lightsON

#this is for circadian counts/index; must define based on how you want to analyze
hrAfter <- 1
hrBefore <- hrAfter

#METADATA READ IN-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#get meta_data, will determine genotypes and loading order
#should have saved "meta_data.csv" with date set up, row, column, genotype, Mutant, N set up, ect.
#read in csv
meta_data <- read.csv("meta_data.csv", header = TRUE, sep=',')
#meta_data$row<-as.integer(meta_data$row) #CHECK IF MONITOR AND ROW ARE LOWERCASE IN META DATA FILE-- OTHERWISE THROWS ERROR

#FFT ANALYSIS AND ANALYSIS_OUTPUT TABLE CONSTRUCTION FROM FFT TABLE--------------------------------------------------------------------------------------------------------------------------------------------------
FFT<-read.csv2("Fourier_Table Data.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)#reads in only the table
#as.integer("string") works, strtoi("string") DOES NOT WORK, returns NA for 08, 09

#function to parse the monitor from the Filename column in Fourier_Table Data.csv. CRUCIAL for analysis_output table
parseM<-function(f)
{
  if (stri_length(f)==11)
  {#if monitor is single digit (1-9)
    m<-as.integer(substr(f, 8, 8))#string values start at 1 not 0
  }
  if (stri_length(f)==12)
  {#if monitor is single digit (10-16)
    m<-as.integer(substr(f, 8, 9))#string values start at 1 not 0
  }
  if (stri_length(f)!=11 && stri_length(f)!=12)
  {
    print("ERROR")
    print(stri_length(f))
  }
  return(m)
}
#function to parse the columns from the Filename column in Fourier_Table Data.csv. CRUCIAL for analysis_output table
parseC<-function(f)
{
  if (stri_length(f)==11)
  {#if monitor is single digit (1-9)
    c<-as.integer(substr(f, 10, 11))
  }
  if (stri_length(f)==12)
  {#if monitor is single digit (10-16)
    c<-as.integer(substr(f, 11, 12))
  }
  if (stri_length(f)!=11 && stri_length(f)!=12)
  {
    print("ERROR")
    print(stri_length(f))
  }
  return(c)
}
#MONITOR, ROW, AND GENOTYPE SHOULD BE ALL LOWER CASE IN META DATA FILE
#parse out the monitor, row, column to include in 
FFT$monitor<-sapply(FFT$Filename, parseM)
FFT$column<-sapply(FFT$Filename, parseC)
FFT$row<-as.integer(((FFT$column-1)/8+1)) #1-8 go to 1, 9-16 go to 2, ect.
c<-c(1,6:8)
#create the Analysis table which will contain all the raw values for each measure for each variable. will be outputted as analysis_output table
Analysis <- FFT[c] #deep copy
Analysis$FFT<-FFT$Amplitude..max.

#PERIODOGRAM ANALYSIS----------------------------------------------------------------------------------------------------------------------------------------------------------------
#PERIOD TIME FRAME LOOKED AT 20-28HR- STANDARD ACCORDING TO SHUBHROZ
#period length (tau) is measured in hours
Periodogram<-read.csv("Periodogram_Table Data.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)#reads in only the table
Analysis$Tau<-Periodogram$Period.1

#PERIODOGRAM RHYTHM ANALYSIS----------------------------------------------------------------------------------------------------------------------------------------------------------------
periodogram_rhythmic<-read.csv2("periodogram_rhythmic.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)#reads in only the table
Analysis$per_rhythmic<-periodogram_rhythmic$rhythmic

#function to define the rhythmicity score of individual fly based on whether I scored it as rhythmic (1) or arrhythmic (0) and its FFT cutoff
RhythmicDefine<-function(f,p)
{
  #define rhythmicity FFT value cutoffs
  #control values quartiled to define out low, medium and high rhythmicity
  #cutoffs determined from our LD data, using known arrhythmic vs rhythimc lines. Not final cuttoffs.
  #0%     25%     50%     75%    100% 
  #.00300 0.00420 0.00610 0.00855 0.01100 
  low<-0.00300
  mid<-0.00420
  high<-0.00610
  if(p==0)
    c<-0
  else if(p==1)
  {
    c<-1
    
    if (f<=low)
    {c<-0
    #all values have to be mathematical integers. 0=arrhythmic
    }
    else if ((f>=low) && (f<mid))
    {c<-1
    #all values have to be mathematical integers. 1=weakly rhythmic
    }
    else if (f>=mid && f<high)
    {c<-2
    #all values have to be mathematical integers. 2=medium rhythmic
    }
    else if (f>=high)
    {c<-3
    #all values have to be mathematical integers. 3= strong rhythmic
    }
    else
      print("ERROR in FFT strength")
  }
  else
    print("ERROR")
  return(c)
}
#define rhythm score
Analysis$RhythmScore<-mapply(RhythmicDefine, Analysis$FFT, p=Analysis$per_rhythmic)

#NUMBER OF FLIES AND NDAY PER LINE ADDED TO META_DATA
#added back to this version (from v19) where the NFlies is part of rhythm score analysis so NFlies will still be calculated with or without rhythm score analysis
for (i in 1:nrow(meta_data))
{
  n=0
  for (j in 1:nrow(Analysis))
  {
    if (meta_data$monitor[i]==Analysis$monitor[j] && meta_data$row[i]==Analysis$row[j]) #THROWS ERROR IF WRONG CASE FOR META DATA ROW/COLUMN- must be lower case
    {
      n=n+1
    }
  }
  meta_data$NDays[i]<-NDays
  meta_data$NFlies[i]<-n
}

#sum each group's rhythmicity by row
for (i in 1:nrow(meta_data))
{
  arrhythmic<-0
  weak<-0
  medium<-0
  strong<-0
  for (j in 1:nrow(Analysis))
  {
    if (meta_data$monitor[i]==Analysis$monitor[j] && meta_data$row[i]==Analysis$row[j]) #THROWS ERROR IF WRONG CASE FOR META DATA ROW/COLUMN- must be lower case
    {
      if (Analysis$RhythmScore[j]==0)
      {
        arrhythmic<- arrhythmic+1
      }
      else if (Analysis$RhythmScore[j]==1)
      { 
        weak<- weak+1
        #all values have to mathematical integers. 1=weakly rhythmic
      }
      else if (Analysis$RhythmScore[j]==2)
      {
        medium<- medium+1
        #all values have to mathematical integers. 2=medium rhythmic
      }
      else if (Analysis$RhythmScore[j]==3)
      {
        strong<- strong+1
        #all values have to mathematical integers. 3= strong rhythmic
      }
      else
        print("ERROR")
    }
  }
  meta_data$Arrhythmic[i]<-arrhythmic
  meta_data$Net_Rhythmic[i]<-(weak+medium+strong)
  meta_data$Weak_Rhythmic[i]<-weak
  meta_data$Mid_Rhythmic[i]<-medium
  meta_data$Strong_Rhythmic[i]<-strong
} 
meta_data$NFlies<-(meta_data$Arrhythmic+meta_data$Net_Rhythmic)

# #COSINOR and NPCRA ACTIVITY PROFILE ANALYSIS------------------------------------------------------------------------------------------------------------------------------------------------------
# activity_profile <- read.csv2("Activity Profile_Table Data.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
# 
# #COSINOR  ANALYSIS
# Analysis$amplitude_COSINOR<-activity_profile$Amplitude
# Analysis$acrophase_COSINOR<-activity_profile$Phase
# Analysis$MESOR_COSINOR<-activity_profile$MESOR
# Analysis$MESOR_Variance_COSINOR<-activity_profile$X..Variance
# 
# #NPCRA
# #onset of M10 and L5
# Analysis$M10_Onset_NPCRA<-activity_profile$M.Start
# Analysis$L5_Onset_NPCRA<-activity_profile$L.Start
# 
# #average amplitude for M10 and L5 and relative amplitude difference between
# Analysis$M10_av_NPCRA<-activity_profile$M.Avg
# Analysis$L5_av_NPCRA<-activity_profile$L.Avg
# Analysis$RA_NPCRA<-activity_profile$RA
# 
# #interdaily stability and intradaily variability
# Analysis$IS_NPCRA<-activity_profile$IS
# Analysis$IV_NPCRA<-activity_profile$IV

# #PHASE SHIFT ANALYSIS------------------------------------------------------------------------------------------------------------------------------------------------------
# #you have to select which of these analyses you to want export and analyze
# #otherwise you'll just have NA data for those columns
# #comment out whichever analyses you don't want
# #this analysis requires both LD and DD data, so this analysis is only done once; as such, it's not included in LD_DD_comp.
# acrophase <- read.csv2("Actogram_Table Data_acrophase.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
# bathyphase <- read.csv2("Actogram_Table Data_bathyphase.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
# onsets <- read.csv2("Actogram_Table Data_onsets.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
# offsets <- read.csv2("Actogram_Table Data_offsets.csv", header = TRUE, sep=',', dec=".", as.is=TRUE, row.names= NULL, stringsAsFactors = FALSE)
# 
# #TAU aka period length
# #each LD tau to Analysis_output
# Analysis$Tau_LD_acro<-acrophase$Tau
# Analysis$Tau_LD_bath<-bathyphase$Tau
# Analysis$Tau_LD_onset<-onsets$Tau
# Analysis$Tau_LD_offset<-offsets$Tau
# 
# #each DD tau to Analysis_output
# Analysis$Tau_DD_acro<-acrophase$Tau.1
# Analysis$Tau_DD_bath<-bathyphase$Tau.1
# Analysis$Tau_DD_onset<-onsets$Tau.1
# Analysis$Tau_DD_offset<-offsets$Tau.1
# 
# #Mean Time of point
# #each LD tau to Analysis_output
# Analysis$Time_LD_acro<-acrophase$Mean.T
# Analysis$Time_LD_bath<-bathyphase$Mean.T
# Analysis$Time_LD_onset<-onsets$Mean.T
# Analysis$Time_LD_offset<-offsets$Mean.T
# 
# #each DD tau to Analysis_output
# Analysis$Time_DD_acro<-acrophase$Mean.T.1
# Analysis$Time_DD_bath<-bathyphase$Mean.T.1
# Analysis$Time_DD_onset<-onsets$Mean.T.1
# Analysis$Time_DD_offset<-offsets$Mean.T.1

#TOTAL COUNTS----DINURAL/NOCTURNAL INDEX-------DAWN AND DUSK ACTIVITY------------------------------------------------------------------------------------------------------------------------------------------------------
countFileNames <- list.files("counts", pattern="*.csv", full.names=TRUE) #lists of all counts files you need to read in to analyze counts- used on loops
#indexing values that get used in counts analysis loop
nightCounts<-0  
dayCounts<-0
Analysis$totalMinSleepPerDay<-0
Analysis$daySleepAv<-0
Analysis$nightSleepAv<-0
Analysis$totalCounts<-0

#column names for the average sleep curves
ns<-list("Filename", "monitor","column","row","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23")
avSleepPerHour<-data.frame(matrix("", ncol=length(ns), nrow = nrow(Analysis)),row.names= rownames(Analysis), stringsAsFactors = FALSE)
colnames(avSleepPerHour)<-ns
#column names for the average counts curves
nc<-list("Filename","monitor","column","row", "0","0.5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10","10.5","11","11.5", "12","12.5","13","13.5","14","14.5","15","15.5","16","16.5","17","17.5","18","18.5","19","19.5","20","20.5","21","21.5","22","22.5","23","23.5")
avCountsPerHalfHour<-data.frame(matrix("", ncol = length(nc), nrow = nrow(Analysis)),row.names= rownames(Analysis))
colnames(avCountsPerHalfHour)<-nc

#building tables for sleep and counts curves
avSleepPerHour[1:4]<-Analysis[1:4]
avCountsPerHalfHour[1:4]<-Analysis[1:4]
avSleepPerHour[5:length(ns)]<-0
avCountsPerHalfHour[5:length(nc)]<-0
sleepBinnedNormalized<-avSleepPerHour
countsBinnedNormalized<-avCountsPerHalfHour


#nested loop for analyzing each count file that was exported from clocklab
#loop iterates over all of the individual counts files
for (i in 1:nrow(Analysis))
{
  #indexing values that get used in counts analysis loop
  nightCounts<-0 #for totaling diurnal and nocturnal index and Light vs Dark activity
  dayCounts<-0
  dawnCounts<-0 #for totaling activity peaks at beginning of day and end of day
  duskCounts<-0
  cntPerMin<-0
  #read in the individual counts file based on the countsFileNames list
  counts<-read.csv(countFileNames[i], header = TRUE, sep=',', dec=".", row.names= NULL, skip=3, stringsAsFactors = FALSE)
  #
  #
  #
  #
  #Changes the name of column 3 in the data frame to Deg so the script will work if the counts per minute data is named deg or cnts/min in the counts file
  colnames(counts)[4] <- "Deg"
  #reads in only the table, no meta data at top of file
  #read.csv reads the csv file in as a table of numerics
  
  #add columns to counts table where assign either 1 or 0 depending on the measure/assay- columns get summed at the end for min of whichever measure
  counts$asleep[1]<-0 #set ahead of time to 0 as won't be iterated over when calculating asleep vs. awake
  counts$asleep[2]<-0 #for determining whether alseep
  counts$asleep[nrow(counts)-1]<-0 #for determining whether alseep
  counts$asleep[nrow(counts)]<-0 #for determining whether alseep
  counts$asleep<-0 #for determining whether alseep
  counts$frag<-0
  counts$boutWake<-0
  counts$boutFallAsleep<-0
  
  #loop iterates over all of the rows of data in the counts file that has been read into counts table to asign asleep vs awake
  for(j in 1:nrow(counts))
  {
    #calculate whether awake or asleep
    if(j>3 && j<(nrow(counts)-2))
    {
      #fly asleep if it's inactive for 5 min
      #get data for 5 min block (2  before and 2 after min of question)
      cntPerMin<-as.numeric(counts$Deg[j])
      cntmin1<-as.numeric(counts$Deg[j-1])
      cntmin2<-as.numeric(counts$Deg[j-2])
      cntplus1<-as.numeric(counts$Deg[j+1])
      cntplus2<-as.numeric(counts$Deg[j+2])
      #add 1 to all points if all the min have no activity (cnt=0)
      if(cntPerMin==0 && cntmin1==0 && cntmin2==0 && cntplus1==0 && cntplus2==0)
      {
        counts$asleep[j]<-1
        counts$asleep[j-1]<-1
        counts$asleep [j-2]<-1
        counts$asleep[j+1]<-1
        counts$asleep[j+2]<-1
      }
      
    }
  }
  
  #loop iterates over all of the rows of data in the counts file that has been read into counts table to determine if fly is experiencing sleep fragmentation
  #dependent upon successful asignment of awake vs asleep (reads $asleep column)
  for(j in 1:nrow(counts))
  {
    if(j>1 && j<(nrow(counts)-1))
    {
      #one definition of sleep fragmentation: sum all times fly woke up for one min then fell asleep again: most accepted measure
      #used to measure sleep bins as well
      
      #calculate when wake for 1 min then fall asleep again immediately
      if(counts$asleep[j]==0 && counts$asleep[j-1]==1 && counts$asleep[j+1]==1)
      {
        counts$frag[j]<-1
      }
      #calculate when wake up
      #DOES NOT remove the sleep fragmentation wakes
      if(counts$asleep[j]==0 && counts$asleep[j-1]==1)
      {
        counts$boutWake[j]<-1
      }
      #calculate when fall asleep
      #DOES NOT remove the sleep fragmentation falling asleep
      if(counts$asleep[j]==0 && counts$asleep[j+1]==1 )
      {
        counts$boutFallAsleep[j]<-(1)
      }
    }
  }   
  
  #bin sleep data into averaged min sleep over Ndays
  x<-0
  cntsum<-sum(counts$Deg, na.rm=TRUE)/NDays #counts/day
  sleepsum<-sum(counts$asleep)/NDays #min sleep/day
  
  #loop calculates sleep per hour and counts per 30 min and creates tables for activity and sleep traces (averaged over Ndays)
  for(k in 0:23)
  {
    #convert EST to zeitgeber time
    y<-k+lightsON
    if (y>(23))
      y<-y-24
    
    cntHrSort<-filter(counts, counts$Hr==y)   #sort by ZT hour
    n<-sum(cntHrSort$asleep) #sum min asleep/hr
    avSleepPerHour[i,k+5]<-as.integer(n)/NDays #average over 5 days and write into table
    sleepBinnedNormalized[i,k+5]<-as.integer(n)/NDays/sleepsum #divide over total sleep and average over 5 days and write into table
    cnt30minSort_First<-filter(cntHrSort, between(cntHrSort$Min,0,29)) #sort data into first half hour
    cnt30minSort_Second<-filter(cntHrSort, between(cntHrSort$Min,30,59)) #sort data into second half hour
    n<-sum(cnt30minSort_First$Deg, na.rm=TRUE)/NDays #sum counts/30 min- first half hour
    avCountsPerHalfHour[i,x+5]<-n #write into table
    countsBinnedNormalized[i,x+5]<-(n/cntsum) #normalized over total activity
    n<-sum(cnt30minSort_Second$Deg, na.rm=TRUE)/NDays #sum counts/30 min- second half hour
    avCountsPerHalfHour[i,x+6]<-n #write into table
    countsBinnedNormalized[i,x+6]<-(n/cntsum) #normalized over total activity
    x<-x+2 #iterate for half hour divisions
  }
  
  #sleep fragmentation- second measure
  #HAMMING INDEX
  x<-1
  hammingDays<-vector(mode="numeric", length=NDays)
  #loop iterates over data and analyzes sleep fragmentation via hamming index
  #iterates over each day- calculates hamming distance for each day (into hammingDays vector)
  #reference?
  for(j in 1:NDays)
  {
    x<-1+1440*(j-1) 
    n<-1440*j
    cntHrSort<-slice(counts,x:n) #subsets data into Zeitberger day 1 (assuming data starts at Zeitberger 0)
    minSleep<-sum(cntHrSort$asleep, na.rm=TRUE)
    n<-(1440-minSleep) #how many potential positions used for hamming index calculation
    cntHrSort$HamIndex<-0
    x<-1
    Hamming<-FALSE
    
    #loop calculates hamming distance
    for (k in 1:n)
    {
      cntHrSort$HamIndex<-0
      cntHrSort$HamIndex[k:(k+minSleep-1)]<-1
      y<-hamming.distance(cntHrSort$asleep, cntHrSort$HamIndex)
      if(!Hamming)
        Hamming<-y
      else if(y<Hamming)
        Hamming<-y
    }
    hammingDays[j]<-Hamming
    
  }
  
  n<-sum(hammingDays, na.rm=TRUE)/NDays/(nrow(counts)-1)
  
  #write measured values into analysis table
  Analysis$HammingDistance[i]<-n 
  Analysis$totalCounts[i]<-sum(counts$Deg, na.rm=TRUE)/NDays
  Analysis$totalMinSleepPerDay[i]<-sum(counts$asleep, na.rm=TRUE)/NDays
  Analysis$NBinPerDay[i]<-(sum(counts$boutWake, na.rm=TRUE)+(sum(counts$boutFallAsleep, na.rm=TRUE)))/NDays/2 #average N fall asleep and N wake incase it falls asleep or wakes at the end/beginning of the data
  #sleep fragmentation is included in NBinsPerDay
  Analysis$FragWake[i]<-sum(counts$frag, na.rm=TRUE)/NDays #N 1 min wake bins
  Analysis$NBinsWithoutFrag[i]<-Analysis$NBinPerDay[i]-Analysis$FragWake[i] #N wake bins without 1 min wakes
  
  #daytime only measures
  #THIS ASSUMES THAT YOUR LIGHTS ON TIME WILL BE BETWEEN 0:00-11:59 EST AND YOUR LIGHTS OFF WILL BE BETWEEN 1:00-23:59 EST 
  cntHrSort<-filter(counts, between(counts$Hr,lightsON,(lightsOFF-1))) #sort for daytime
  
  Analysis$dayCounts[i]<-sum(cntHrSort$Deg, na.rm=TRUE)/NDays #calculates average total counts during light/day
  Analysis$daySleepAv[i]<-sum(cntHrSort$asleep, na.rm=TRUE)/NDays #calculates average total sleep during light/day
  Analysis$dayBins[i]<-(sum(cntHrSort$boutWake, na.rm=TRUE)+(sum(cntHrSort$boutFallAsleep, na.rm=TRUE)))/NDays/2 #N daytime sleep bins
  Analysis$dayFrag[i]<-sum(cntHrSort$frag, na.rm=TRUE)/NDays #daytime sleep fragmentation
  Analysis$dayBinsWithoutFrag[i]<-Analysis$dayBins[i]-Analysis$dayFrag[i] #daytime N wake bins without 1 min wakes
  
  #nighttime only measures
  #THIS ASSUMES THAT YOUR LIGHTS ON TIME WILL BE BETWEEN 0:00-11:59 EST AND YOUR LIGHTS OFF WILL BE BETWEEN 12:00-23:59 EST
  cntHrSort<-filter(counts, between(counts$Hr,lightsOFF, 24)) #sort for lights off to midnight
  cntHrSort2<-filter(counts, between(counts$Hr,0, (lightsON-1))) #sort for midnight to lights on
  cntHrSort3<-rbind(cntHrSort, cntHrSort2) #combine above to get for lights off periodf
  
  Analysis$nightCounts[i]<-(sum(cntHrSort3$Deg, na.rm=TRUE))/NDays #calculates average total counts during dark
  Analysis$nightSleepAv[i]<-(sum(cntHrSort3$asleep, na.rm=TRUE))/NDays #sleep during dark
  Analysis$nightBins[i]<-(sum(cntHrSort3$boutWake, na.rm=TRUE)+(sum(cntHrSort3$boutFallAsleep, na.rm=TRUE)))/NDays/2 #sleep bins at night
  Analysis$nightFrag[i]<-sum(cntHrSort3$frag, na.rm=TRUE)/NDays #sleep frag at night
  Analysis$nightBinsWithoutFrag[i]<-Analysis$nightBins[i]-Analysis$nightFrag[i] #daytime N wake bins without 1 min wakes
  
  #sort for times when activity peaks during light switching times
  cntHrSort<-filter(counts, between(counts$Hr,(lightsOFF-hrBefore), (lightsOFF+hrAfter-1)))
  Analysis$duskCounts[i]<-sum(cntHrSort$Deg, na.rm=TRUE)/NDays#calculates average total counts N hours before and N hours after lights off/day
  cntHrSort<-filter(counts, between(counts$Hr,(lightsON-hrBefore), (lightsON+hrAfter-1)))
  Analysis$dawnCounts[i]<-sum(cntHrSort$Deg, na.rm=TRUE)/NDays #calculates average total counts N hours before and N hours after lights on/day
} #end nested loops

#add remaining measures that don't need to be iterated over (are already calculable from values in analysis table) to analysis table
Analysis$UrnalIndex_Counts<- (Analysis$dayCounts-Analysis$nightCounts)/Analysis$totalCounts
Analysis$circadianIndex_Counts<- (Analysis$dawnCounts + Analysis$duskCounts)/Analysis$totalCounts
Analysis$BinLength_Min<-Analysis$totalMinSleepPerDay/Analysis$NBinPerDay
Analysis$BinLengthWithoutFrag<-(Analysis$totalMinSleepPerDay+Analysis$FragWake)/Analysis$NBinsWithoutFrag
Analysis$dayBinLength<-Analysis$daySleepAv/Analysis$dayBins
Analysis$dayBinLengthWithoutFrag<-(Analysis$daySleepAv+Analysis$dayFrag)/Analysis$dayBinsWithoutFrag
Analysis$nightBinLength<-Analysis$nightSleepAv/Analysis$nightBins
Analysis$nightBinLengthWithoutFrag<-(Analysis$nightSleepAv+Analysis$nightFrag)/Analysis$nightBinsWithoutFrag
Analysis$dayCounts_Ratio<-Analysis$dayCounts/Analysis$totalCounts
Analysis$nightCounts_Ratio<-Analysis$nightCounts/Analysis$totalCounts
Analysis$daySleep_Ratio<-Analysis$daySleep/Analysis$totalMinSleepPerDay
Analysis$nightSleep_Ratio<-Analysis$nightSleep/Analysis$totalMinSleepPerDay

#WRITE ALL FILES----------------------------------------------------------------------------------------------------------------------------------------------------------------
#create subfolder to put analysis into (if one doesn't already exist)
if(!dir.exists(analysis_folder))
  dir.create(analysis_folder, showWarnings = TRUE, recursive = FALSE)
setwd(analysis_folder)
write.csv(analysis_meta, file="date and version.csv")
write.csv(Analysis, file= "analysis_output.csv") #this works fine just leave it alone, dont use csv2 or will write as ; and excel throws a shit fit and clarifying divider doesnt help
write.csv(meta_data, file= "meta_data_analysis.csv") 
write.csv(avSleepPerHour, file= "average_sleep_hour_bin.csv")
write.csv(avCountsPerHalfHour, file= "average_counts_half_hour_bin.csv")
write.csv(sleepBinnedNormalized, file= "average_sleep_hour_bin_normalized.csv")
write.csv(countsBinnedNormalized, file= "average_counts_half_hour_bin_normalized.csv")
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
