# ##############################################################################
# cBioPortal CNV data analysis
# ##############################################################################
# remove all objects
rm(list = ls(all = TRUE))

# load the required libraries
library(data.table) # to load very large data files
library(survival)
library(ggplot2)
library(ggthemes)
library(survminer) # required for ggsurvplot

# go to working directory 
setwd("~/Dropbox/Manuscripts/Project X/6-CNV/cBioportal CNV analysis/")


# import Excel file
library("openxlsx")
# mRNA copy number adjusted by NF				
CNV.DF <- read.xlsx("cBioPortal workings.xlsx", sheet = 2, startRow = 1, colNames = TRUE)  # CNV data
# sheet 1: sampleinfo
# sheet 2: CNVs
# sheet 3: sample type info. (metastatic or not)
# sheet 4: OS (overall survival)
# sheet 5: DFS (disease-free survival)
head(CNV.DF)
colnames(CNV.DF)
class(CNV.DF) # df
# #####################
# RENAME CNV TYPE
# change HETLOSS; to CNV1
# change HOMDEL; to CNV;
table(CNV.DF$CIC.event)
# HETLOSS;  HOMDEL; 
#  148       29 

CNV.DF$CIC.event <- gsub("HETLOSS;",1,CNV.DF$CIC.event,ignore.case=T)
CNV.DF$CIC.event <- gsub("HOMDEL;",0,CNV.DF$CIC.event,ignore.case=T)
# convert NA to normal gene copy number = 2
CNV.DF[c("CIC.event")][is.na(CNV.DF[c("CIC.event")])] <- 2
CNV.DF$CIC.event
#
#
# 0    1    2 
# 29  148 1136 

table(CNV.DF$dataset)


# #####################

# #####################
# LOAD SAMPLETYPE INFO.
# #####################
sampletype.DF <- read.xlsx("cBioPortal workings.xlsx", sheet = 3, startRow = 1, colNames = TRUE)  # CNV data
head(sampletype.DF)
table(sampletype.DF$Sample.Type)
# Metastasis        Metastatic           Primary           PRIMARY     Primary Tumor Primary/Localized 
# 369                88               202               179               498                19 
# use 
# Primary
# Metastatic
sampletype.DF <- as.data.frame( lapply(sampletype.DF, function(x) {gsub("Metastasis", "Metastatic", x)}) )
sampletype.DF <- as.data.frame( lapply(sampletype.DF, function(x) {gsub("PRIMARY", "Primary", x)}) ) 
sampletype.DF <- as.data.frame( lapply(sampletype.DF, function(x) {gsub("Primary Tumor", "Primary", x)}) )
sampletype.DF <- as.data.frame(  lapply(sampletype.DF, function(x) {gsub("Primary/Localized", "Primary", x)}) )
table(sampletype.DF$Sample.Type)


#   Metastatic    Primary 
#     457        898 

# #####################
# ONLY KEEP SAMPLES IN CNV.DF 
# #####################
class(sampletype.DF) # df 
class(CNV.DF) # df 
colnames(sampletype.DF) # df 
colnames(CNV.DF) # df 
table(sampletype.DF$data.set) # NEPC (Trento/Cornell/Broad 2016) still there
table(CNV.DF$dataset) # NEPC (Trento/Cornell/Broad 2016) still there

mergedDF <- merge(CNV.DF,sampletype.DF, by.x='Case.ID', by.y='Patient.ID',all=FALSE) # all=FALSE ... do not keep
head(mergedDF)

metastatic.DF <- mergedDF[which( mergedDF$Sample.Type == "Metastatic"),]
primary.DF <- mergedDF[which( mergedDF$Sample.Type == "Primary"),]

table(metastatic.DF$CIC.event)
# Met
(10/(10+95+337))*100 # CNV0 2.3%
(95/(10+95+337))*100 # CNV1 21.5%
(337/(10+95+337))*100 # normal 76.2%

# Prim.
table(primary.DF$CIC.event)
(19/(19+53+797))*100 # 2.2%
(53/(19+53+797))*100 # 6.1%
(797/(19+53+797))*100 # 91.7%
# CIC: HETLOSS HOMDEL;


# how many of each type (for data set)
colnames(mergedDF)
table(mergedDF$dataset)
# NEPC (Trento/Cornell/Broad 2016)    Prostate (Broad/Cornell 2013)           Prostate (FHCRC, 2016)                  Prostate (MICH) 
# 107                               56                              149                               61 
# Prostate (MSKCC 2010)            Prostate (MSKCC 2014)                  Prostate (SU2C)                  Prostate (TCGA) 
# 194                              104                              150                              492 
datasetofinterest <- "Prostate (TCGA)"
Prims <- length(which(mergedDF$dataset == datasetofinterest & mergedDF$Sample.Type == "Primary"))  # Primary 34
Mets <- length(which(mergedDF$dataset == datasetofinterest & mergedDF$Sample.Type == "Metastatic"))  # Metastatic 73
paste("(", Prims," primary and ",Mets," metastatic tumors)", sep="")





# graphing

# METASTATIC
WTcount.WT <-   length(which(metastatic.DF$CIC.event == "2")) # 185
WTcount.HETLOSS <-   length(which(metastatic.DF$CIC.event == "1")) # 50
WTcount.HOMDEL <-   length(which(metastatic.DF$CIC.event == "0")) # 4
df <- data.frame(
  group = c("WT", "HETLOSS", "HOMDEL"),
  value = c(WTcount.WT, WTcount.HETLOSS, WTcount.HOMDEL)
)
head(df)
library(ggplot2)
# Barplot
bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0)
pie
# Use custom color palettes
pie + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme_void()

pie + scale_fill_manual(values=c("#6f4295","#ce3629","#e0922e")) + theme_void() # Nature MS
# orange e0922e   red  ce3629

# WT HETLOSS HOMDEL
pdf(paste("CIC-alterations-metastatic samples.pdf",sep=""))
pie + scale_fill_manual(values=c("#61b983","#6774fa","#e0922e")) + theme_void()
print(p)
dev.off()
# LNCaP green 
# 61b983
# PC3 blue 
# 6774fa

# PRIMARY
WTcount.WT <-   length(which(primary.DF$CIC.event == "2")) # 185
WTcount.HETLOSS <-   length(which(primary.DF$CIC.event == "1")) # 50
WTcount.HOMDEL <-   length(which(primary.DF$CIC.event == "0")) # 4
df <- data.frame(
  group = c("WT", "HETLOSS", "HOMDEL"),
  value = c(WTcount.WT, WTcount.HETLOSS, WTcount.HOMDEL)
)
head(df)
library(ggplot2)
# Barplot
bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0)
pie
# Use custom color palettes
pie + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme_void()

pie + scale_fill_manual(values=c("#6f4295","#ce3629","#e0922e")) + theme_void() # Nature MS
# orange e0922e   red  ce3629

# WT HETLOSS HOMDEL
pdf(paste("CIC-alterations-primaries samples.pdf",sep=""))
pie + scale_fill_manual(values=c("#61b983","#6774fa","#e0922e")) + theme_void()
print(p)
dev.off()





# ##############################################################################################################
# ADD SURVIVAL INFORMATION
# ##############################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOT MANY DATA SETS WITH SURVIVAL INFORMATION!


# relapse info.
dfs.DF <- read.xlsx("cBioPortal workings.xlsx", sheet = 5, startRow = 1, colNames = TRUE)  # CNV data
head(dfs.DF)
names(dfs.DF)[names(dfs.DF) == 'Disease-free.Survival.(Months)'] <- 'DFS'
names(dfs.DF)[names(dfs.DF) == 'Disease-free.Survival.Status'] <- 'DFS.EVENT'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# DO NOT USE THE TCGA DFS/OS INFO. IN THE TCGA cBioPortal
# This script will load *updated* TCGA clinical information from cBioPortal (continuously-updated)
dataset <- 'PRAD'
source("2D-clinical information.R")
head(clinicalDF)
#  Sample.ID   Sample.Type   DFS   DFS.EVENT
clinicalDF$Sample.ID <- gsub("-01","",clinicalDF$Sample.ID,ignore.case=T)
names(clinicalDF)[names(clinicalDF) == 'Sample.ID'] <- 'Patient.ID'
#
head(clinicalDF)
#  Sample.ID   Sample.Type   DFS   DFS.EVENT
clinicalDF <- subset(clinicalDF, select=c(Patient.ID,DFS,DFS.EVENT))
# just rbind them since they have the same order
# colnames 
clinicalDF$data.set <- "Prostate (TCGA)"
colnames(clinicalDF)
colnames(dfs.DF)

dfs.DF.final <- rbind(dfs.DF,clinicalDF)

table(dfs.DF.final$DFS.EVENT)
#  DiseaseFree            Recurred Recurred/Progressed 
#    616                  84                  92 
dfs.DF.final[dfs.DF.final=="DiseaseFree"]<-0
dfs.DF.final[dfs.DF.final=="Recurred/Progressed"]<-1
dfs.DF.final[dfs.DF.final=="Recurred"]<-1
dfs.DF.final$DFS <- as.numeric(dfs.DF.final$DFS)
dfs.DF.final$DFS.EVENT <- as.numeric(dfs.DF.final$DFS.EVENT)
table(dfs.DF.final$DFS.EVENT)
#


mergedDF.dfs <- merge(mergedDF,dfs.DF.final, by.x='Case.ID', by.y='Patient.ID',all=FALSE) # all=FALSE ... do not keep
head(mergedDF.dfs)
mergedDF.dfs$CIC.event <- as.numeric(mergedDF.dfs$CIC.event) 
class(mergedDF.dfs$CIC.event)

# colnames(mergedDF.dfs)
mergedDF.dfs <- subset(mergedDF.dfs, select=c(Case.ID,DFS,DFS.EVENT,CIC.event,Sample.Type,dataset))

# if we want to only examine metastatic samples
length(which(tempdata$Sample.Type == "Metastatic")) # 41  
length(which(tempdata$Sample.Type == "Primary"))  # 750
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ##############################################################################################################
# SUBSET SAMPLES HERE
# ##############################################################################################################
# 0   1   2 
# 18  46 727 
#  0   1   2 
# 18  46 727 
which(mergedDF.dfs$CIC.event == 0 & mergedDF.dfs$Sample.Type == "Metastatic") # only one. sample 186
which(mergedDF.dfs$CIC.event == 1 & mergedDF.dfs$Sample.Type == "Metastatic") # 6

# try various combinations here
analysisname <- "Metastatic-only-WTvsHETvsHOM.pdf"
tempdata <- mergedDF.dfs[which(mergedDF.dfs$Sample.Type == "Metastatic"),]
tempdata <- mergedDF.dfs[which(mergedDF.dfs$Sample.Type != "Metastatic"),] # PRIMARY KEPT
tempdata <- tempdata[which(tempdata$CIC.event != 1),] # only keep HOMDEL
tempdata <- tempdata[which(tempdata$CIC.event != 0),] # only keep HETLOSS
source("SUB_grouping.R")








# ##############################################################
# survival analysis
# ##############################################################

# separate, so we obtain P-values for each condition


# METASTATIC--WT vs HET
tempdata <- mergedDF.dfs
tempdata <- mergedDF.dfs[which(mergedDF.dfs$Sample.Type == "Metastatic"),] # PRIMARY KEPT
tempdata <- tempdata[which(tempdata$CIC.event != 0),] # only keep HETLOSS
try(source("SUB_grouping.R"))
# METASTATIC--WT vs HOM
tempdata <- mergedDF.dfs
tempdata <- mergedDF.dfs[which(mergedDF.dfs$Sample.Type == "Metastatic"),] # PRIMARY KEPT
tempdata <- tempdata[which(tempdata$CIC.event != 1),] # only keep HOMDEL
try(source("SUB_grouping.R"))

# PRIMARY--WT vs HET
tempdata <- mergedDF.dfs
tempdata <- mergedDF.dfs[which(mergedDF.dfs$Sample.Type != "Metastatic"),] # PRIMARY KEPT
tempdata <- tempdata[which(tempdata$CIC.event != 0),] # only keep HETLOSS
try(source("SUB_grouping.R"))
# PRIMARY--WT vs HOM
tempdata <- mergedDF.dfs
tempdata <- mergedDF.dfs[which(mergedDF.dfs$Sample.Type != "Metastatic"),] # PRIMARY KEPT
tempdata <- tempdata[which(tempdata$CIC.event != 1),] # only keep HOMDEL
try(source("SUB_grouping.R"))

survdiff(Surv(DFS, DFS.EVENT) ~ group ,data=dummyDF , rho=0)  
#
DFSsurv <- survdiff(Surv(DFS, DFS.EVENT) ~ group ,data=dummyDF, rho=0)
p.val.DFS <- 1 - pchisq(DFSsurv$chisq, length(DFSsurv$n) - 1)
signaturePvalue <- p.val.DFS



# for figure with all three
# normal vs heterozygous deletion (and include HOM)
analysisname <- "WTvsHETvsHOM-Primary-v2.pdf"
tempdata <- mergedDF.dfs
source("SUB_grouping.R")

#
#
# can draw several groups here
oldpalette <- c("#6f4295","#ce3629","#e0922e")
currentpalette <- c("#61b983","#6774fa","#e0922e")
fit=survfit(Surv(DFS, DFS.EVENT) ~ group,data=dummyDF)
pdf(paste(analysisname,"-survival-DFS.pdf",sep=""))
p <- ggsurvplot(fit, risk.table = FALSE,
                pval = TRUE, risk.table.y.text.col = TRUE, break.time.by=1,conf.int=FALSE, censor=TRUE, legend="top", palette=currentpalette)
print(p)
dev.off()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
