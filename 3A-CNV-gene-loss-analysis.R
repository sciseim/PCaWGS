## 3A-CNV-gene-loss-analysis.R

# Load the cn.mops package
library(cn.mops)


# run on Ewha computer since I do not have enough space and RAM on my Mac
#
# before runnign, touch all bam and bai files to ensure no problems with bam index time stamps
#@ touch /media/sciseim/mogwai/PCaWGS_MS/CNV/BAM/*.ba*

# BAM files (and BAI)
# /media/sciseim/mogwai/PCaWGS_MS/CNV/BAM
# PC3-reads.bam 
# LNCaP-reads.bam
# NP-reads.bam
# all >100GB and aligned to hg38
# *PAIRED* READS

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Please denote which of the samples are interested in
sampleofinterest <- "LNCaP"
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


## set path/to/download/directory, e.g.,
destdir <- "/media/sciseim/mogwai/PCaWGS_MS/CNV/BAM"
setwd(destdir)

bamFiles <- file.path(destdir,
                      c(paste(sampleofinterest,"-reads.bam",sep=""), "NP-reads.bam"))

# And for what it's worth the sequence names ('levels') in the bam files can be discovered with 
library(Rsamtools)
seqlevels(BamFileList(bamFiles))


# refSeqName=c(paste0("",1:22),"X","Y")
## 2. We can bin and count the reads
# no chr prexix for hg38 --e.g. Y, not chrY
# The parameter "refSeqName"  relates to the names of the chromosomes as they are named in the BAM file header. Typically, # these are "chr1", "chr2",... or "1", "2",... depending on the version of the reference genome.
reads_gr <- getReadCountsFromBAM(bamFiles, sampleNames = c("tumor", "normal"),refSeqName=c(paste0("",1:22),"X","Y"), WL = 10000, mode = "paired")  
# test using chrY (tiny)
# worked

 

## 3. Normalization
## We need a special normalization because the tumor has many large CNVs
# 
# Normalize quantitative NGS data in order to make counts comparable over samples. Scales each samples’ reads such that the # coverage is even for all samples after normalization.
X <- normalizeGenome(reads_gr, normType="poisson")

# save output
saveRDS(X,file=paste("X.normalized-",sampleofinterest,"vsNP.Rds",sep="")) # normalised

# CN2 is the normal copy number for diploid samples. CN1 is a heterozygous deletion and CN0 is a homozygous deletion
#

## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("cn.mops")
library("cn.mops")
library(GenomeInfoDb) # should be autoloaded with cn.mops
library(Rsamtools)

# library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# library("org.Hs.eg.db")
library("DNAcopy")
library("biomaRt")
## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("DNAcopy")
setwd("~/Dropbox/Manuscripts/Project X/6-CNV/cn.mops")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Please denote which of the samples are interested in
sampleofinterest <- "PC3"
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# load CNV data
X <- readRDS(paste("X.normalized-",sampleofinterest,"vsNP.Rds",sep=""))


## 4. Detect cnv's
# here 1: tumor 2: normal in X
# lowerThreshold = -0.9 is default
# priorImpact = 1 is default
resRef <- referencecn.mops(cases=X[,1],controls=X[,2],classes=c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6","CN7","CN8","CN16","CN32","CN64","CN128"),I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64),segAlgorithm="DNAcopy",lowerThreshold=  -0.9, priorImpact = 1 )
resRef <- calcIntegerCopyNumbers(resRef)
(cnvs(resRef))
# MY OLD CODE
# ref_analysis <- referencecn.mops(X[,1], X[,2],
#                                  norm=0, 
#                                  I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64), 
#                                  classes = paste0("CN", c(0:8, 16, 32, 64, 128)),
#                                  segAlgorithm="DNAcopy")
resRef <- calcIntegerCopyNumbers(resRef)
(cnvs(resRef))
ref_analysis <- resRef



# save/load
 saveRDS(ref_analysis,paste("resRef-",sampleofinterest,"vsNP.Rds",sep=""))
 ref_analysis <- readRDS(paste("resRef-",sampleofinterest,"vsNP.Rds",sep=""))



## Analyzing: Sample.1
resCNMOPS <- calcIntegerCopyNumbers(ref_analysis)

## 5. Visualize the cnv's
pdf(paste("visualise-CNVs-",sampleofinterest,".PDF",sep=""),width=50,height=50)
h <- segplot(resCNMOPS)
print(h)
dev.off()


# Here the x-axis represents the genomic position and the y-axis represents the log ratio of read counts and copy number 
# call of each segment (red)
remove(human_cn)
human_cn <- cnvr(resCNMOPS)
human_cn

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# NEXT --FIX THE SCRIPT --OBTAIN ALL GENES IN A REGION OF INTEREST!
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# CDH18 	 cadherin 18 [Source:HGNC Symbol;Acc:HGNC:1757] 	 CN0 	 chr5 	 - 	19850001	19960000	110000	 Q13634
# 110,000 bp in chr5:19850001-19960000 THIS IS JUST THE GENE COORDINATES, BUT THE WIDTH SAY=25,000 bp


TESThuman_cnDF <- as.data.frame(human_cn)
head(TESThuman_cnDF)
tail(TESThuman_cnDF)
# 312 

dat <- TESThuman_cnDF
dat$temp <- do.call(paste, c(dat[c("start","end")], sep = ":")) 
names(dat)[names(dat)=="seqnames"] <- "chr"
names(dat)[names(dat)=="tumor"] <- "CNV"
#
dat$filterlist <- do.call(paste, c(dat[c("chr", "temp")], sep = ":")) 
# dat$filterlist <- do.call(paste, c(dat[c("chr", "start","end")], sep = ":")) 
dat$filterlist <- gsub('chr', '', dat$filterlist) # need to remove chr again for it to work with this genome


dat <- dat[which(dat$CNV == "CN0"),] # keep homozygous deletions only 
# PC3 = 21
try(source("FUN.biomaRt.CNV.regions.R"))

# CN1  
# dat <- dat[which(dat$CNV == "CN1"),] # keep heterozygous deletions only
# CN1 is a heterozygous deletion


# tail(dat)
# e.g.    chr    start      end   width strand CNV              temp          filterlist
# 510   Y  9450001 10270000  820000      * CN0  9450001:10270000  Y:9450001:10270000

# now have
biomaRtLookUp.subset # genes from biomaRt *with* 'DATROW' (5:19850001:19960000) that we can use to match up with dat 'filterlist'
dat

mergedDF <- merge(dat,biomaRtLookUp.subset, by.x=c("filterlist"),by.y=c("DATROW"),all = FALSE)
head(mergedDF)                 

# change order
mergedDF <- mergedDF[,c("hgnc_symbol","description", "CNV", "chr", "filterlist","start_position", "end_position" )]
# re-name some columns
names(mergedDF)[names(mergedDF)=="hgnc_symbol"] <- "symbol"
names(mergedDF)[names(mergedDF)=="filterlist"] <- "CNV region"
names(mergedDF)[names(mergedDF)=="start_position"] <- "gene start"
names(mergedDF)[names(mergedDF)=="end_position"] <- "gene end"
head(mergedDF)
# get rid of [Source etc]  ... lipase F, gastric type [Source:HGNC Symbol;Acc:HGNC:6622] 
mergedDF <- mergedDF
mergedDF$description <-   sub("\\[.*","",mergedDF$description)
head(mergedDF)
# write output
write.table(mergedDF, file = paste(sampleofinterest,"-hg38-CNV.txt",sep=""), sep = " \t ", row.names = FALSE, quote = FALSE)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# next, let us look at RNAseq 


try(source("3B-load.RNAseq.R"))
try(source("3C-merge.RNAseq.and.CNV.R"))


head(mergedDF.final.UPCfiltered)
# still does not include non-coding RNAs. Need to use featureCounts with '-t exon' , rather '-t CDS'
# i.e. count all exons, not just CDS exons
# checking the output counts
# MIR548AT = ENSG00000264314
# grep "ENSG00000264314" *  works



# data in mergedDF.final.UPCfiltered
# ##########################
# WRITE OUT A TABLE 
# ##########################
# paste(sampleofinterest,"putative-lost-genes.txt",sep="_")
# write.table(mergedDF.final, "putative-lost-genes-PC3-assembly.txt", row.names=F, col.names=T, sep="\t", quote=F)


# change the CNV region to a more common syntax (ensembl biomaRt region before)
mergedDF.final.UPCfiltered$`CNV region` <- gsub("(.*)\\:", "\\1-\\2", mergedDF.final.UPCfiltered$`CNV region`)
mergedDF.final$`CNV region` <- gsub("(.*)\\:", "\\1-\\2", mergedDF.final$`CNV region`)



write.table(mergedDF.final.UPCfiltered, paste(sampleofinterest,"putative-lost-genes.txt",sep="_"), row.names=F, col.names=T, sep="\t", quote=F)

# and save one without filtering out genes that are not expressed by NP (e.g. no endogenous expression)

write.table(mergedDF.final, paste(sampleofinterest,"putative-lost-genes-noUPCfiltering.txt",sep="_"), row.names=F, col.names=T, sep="\t", quote=F)