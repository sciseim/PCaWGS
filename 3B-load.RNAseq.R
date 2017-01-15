# 2-RNAseq

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ##########################
# load RNAseq data
# ##########################
matching <- read.delim("./5-RNAseq/samplelist.txt", header=F, comment.char = "#")

result = list()
for (i in 1:nrow(matching)) {
  
  temp <- read.delim(paste("./5-RNAseq/count/", matching[i,1], ".count", sep=""), comment.char = "#", header=T)
  
  # get rid of any columns where gene name is missing
  #@ temp <- temp[!is.na(temp$Geneid),]
  
  # temp[,1] = gsub(".*synthetic_", "", temp[,1])
  rownames(temp) <- temp[,1]
  
  head(temp)
  tail(temp)
  rownames(temp)
  
  temp <- temp[,c(1,7)]
  
  #  temp = temp[,c(1,2,6,7)] to re-order
  #  colnames(temp) = c("Ortholog.ID", "Species.ID", "Length", "Count")
  
  result[[as.character(matching[i,1])]] <- temp
}


# which gene names are common?
gene.common <- table(unlist(lapply(result, rownames)))
nosamples <- length(result)
head(gene.common)
gene.common <- names(gene.common)[gene.common==nosamples]  # change to # of samples!  
head(gene.common)

# only keep common genes in the list results (really a dfList)
result.common <- lapply(result, function(x) x[gene.common,])
# head(result.common)

# create the final table with counts 
count.table <- do.call(cbind.data.frame, lapply(result.common, function(x) x[,2])) # ,2 is count
rownames(count.table) <- gene.common

head(count.table)
tail(count.table)

# NORMALISE AND FILTER
# Normalization time!
library("edgeR")

tmm <- calcNormFactors(count.table, method = "TMM")
lib.size <- colSums(count.table)/10^6
eff.lib.size <- lib.size * tmm
count.table.scale <- t(t(count.table) / eff.lib.size)

min.count = apply(count.table.scale, 1, function(x) sum(x<3))
table(min.count)
head(count.table.scale)

# get rid of rows with only zeros...
 count.table.scale.no0 <- count.table.scale[ rowSums(count.table.scale)!=0, ] 
 head(count.table.scale.no0)
 count.table.scale <- count.table.scale.no0 
 head(count.table.scale)

# replace 0s with NA
count.table.scale.nona <- count.table.scale
count.table.scale.nona[count.table.scale.nona==0] = NA
# count.table.scale.nona[count.table.scale.nona==0] = 0.0000001  # Better! 
# expect many to be 0 prior to hypoxia insult!
count.table.scale <- count.table.scale.nona
head(count.table.scale)

# number of orthologs
length(count.table.scale)/nosamples # 18,744
counts <- count.table.scale
class(counts) # matrix

head(counts)
counts <- as.data.frame(counts)
counts$ensembl.gene <- row.names(counts)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FFAR3 = ENSG00000185897
counts[which(counts$ensembl.gene == "ENSG00000185897"),]

# ##########################
# SCAN.UPC
# ##########################
# load the required libraries
library(SCAN.UPC)

class(counts) # df
# classmatrix <- as.matrix(subset(counts, select=-c(ensembl.gene)))
classmatrix <- (subset(counts, select=-c(ensembl.gene)))
classmatrix[is.na(classmatrix)] <- 0
classmatrix <- as.matrix(classmatrix)

# run UPC
dataset <- "PC3-reads"
UPCnormalisedname <- paste("UPCnormalised.",dataset,sep="")
# UPC_Generic
# Generic function to apply Universal exPression Codes (UPC) trans- formation
# This function can be used to derive UPC values to any type of gene-expression data. It requires the
# user to specify expression values for many genes (or transcripts, exons, or probes). And optionally,
# the user can specify the length and/or GC content (proportion of G or C bases) for the corresponding
# genomic region (e.g., gene). If these values are specified, the UPC algorithm corrects for biases
# resulting from length or GC content.
# 

# MAY REMOVE BELOW
# dim(classmatrix[2])
# # run for each column
# head(classmatrix)  
# colnames(classmatrix)
# for(i in 1:dim(classmatrix)[1]){
#   #  print("YOKI!")
#   print(colnames(classmatrix)[1])
#   UPC.scores <- UPC_Generic(dat$E[,i], verbose = FALSE)
# }
# #  MAY REMOVE ABOVE



head(classmatrix)
remove(UPC.analysis.list)
UPC.analysis.list <- list()
# need to run UPC for each sample
for (i in 1:length((colnames(classmatrix))) )  # 9 samples
{
  # i <- 1 # dummy
  print(colnames(classmatrix)[i])
  samplevalues <- as.matrix((classmatrix)[,i])                 # all the values from column 1
  length(samplevalues) # 18,744
  UPC.scores <- UPC_Generic(samplevalues, verbose = FALSE) 
  UPC.analysis.list[[i]] <- UPC.scores    # save for each i (here: TCGA data set) 
}

UPC.analysis.DF <- do.call(rbind.data.frame, UPC.analysis.list)  # convert the list to a data frame!
class(UPC.analysis.DF) # df
head(UPC.analysis.DF)
UPC.analysis.DF <- as.data.frame(t(UPC.analysis.DF)) # force df, or default to matrix
# 2 18744
row.names(UPC.analysis.DF) <- row.names(counts)
head(UPC.analysis.DF)
class(UPC.analysis.DF) # df
UPC.analysis.DF$ensembl.gene <- row.names(UPC.analysis.DF)
colnames(UPC.analysis.DF) <- colnames(counts)
head(UPC.analysis.DF)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
filterlist <- UPC.analysis.DF$ensembl.gene

# look up gene symbol
# biomaRt look-up
ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

infotime=getBM(attributes = c("hgnc_symbol","description","entrezgene", "chromosome_name", "start_position", "end_position","gene_biotype","ensembl_gene_id"),filters = c("ensembl_gene_id"),values = list(filterlist), mart = ensembl)
#head(infotime)
# hgnc_symbol   
# match by ensembl_gene_id

UPC.analysis.DF.symbolnames <- merge(UPC.analysis.DF,infotime, by.x=c("ensembl.gene"),by.y=c("ensembl_gene_id"),all = FALSE)
# colnames(UPC.analysis.DF.symbolnames)    
UPC.analysis.DF.symbolnames <- subset(UPC.analysis.DF.symbolnames, select=-c(description,entrezgene,chromosome_name,start_position,end_position,gene_biotype, ensembl.gene)) # -c , so removing columns
# head(UPC.analysis.DF.symbolnames)
