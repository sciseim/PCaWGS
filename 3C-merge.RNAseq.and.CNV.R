
# ##########################
# merge mergedDF with RNAseq data
# ##########################
head(mergedDF)

# if no UPC (normalised counts/CPM)
#mergedDF.final <- merge(mergedDF, counts, by.x = "ensembl_gene_id", by.y = "ensembl.gene")

# if prefer to use UPC values
mergedDF.final <- merge(mergedDF, UPC.analysis.DF.symbolnames, by.x = "symbol", by.y = "hgnc_symbol")


head(mergedDF.final)
colnames(mergedDF.final)

# keep RefSeq ID only since CDS and AA referes to these
# mergedDF.final <- mergedDF.final[,c("symbol","description","chr","CNV","CNV region",
#                                     
#                                     
#                                     "End","Consequence","Exon","RefSeq ID","CDS mutation","AA mutation","COSMIC ID","COSMIC occurrence","120215-UNC10-SN254-0327-AC0CMCACXX-ACTTGA-L005", "120502-UNC14-SN744-0235-BD0YUTACXX-ACTTGA-L005", "130221-UNC9-SN296-0338-BC1PYCACXX-TGACCA-L008", "SRR1735560", "SRR1735559", "SRR1735558", "SRR1767372", "SRR1767371", "SRR1767370")]

# give the RNAseq sample sane names
names(mergedDF.final)[names(mergedDF.final) == '120215-UNC10-SN254-0327-AC0CMCACXX-ACTTGA-L005'] <- 'NP.1'
names(mergedDF.final)[names(mergedDF.final) == '120502-UNC14-SN744-0235-BD0YUTACXX-ACTTGA-L005'] <- 'NP.2'
names(mergedDF.final)[names(mergedDF.final) == '130221-UNC9-SN296-0338-BC1PYCACXX-TGACCA-L008'] <- 'NP.3'
#
names(mergedDF.final)[names(mergedDF.final) == 'SRR1735560'] <- 'LNCaP.1'
names(mergedDF.final)[names(mergedDF.final) == 'SRR1735559'] <- 'LNCaP.2'
names(mergedDF.final)[names(mergedDF.final) == 'SRR1735558'] <- 'LNCaP.3'
#
names(mergedDF.final)[names(mergedDF.final) == 'SRR1767371'] <- 'PC3.1'
names(mergedDF.final)[names(mergedDF.final) == 'SRR1767372'] <- 'PC3.2'
names(mergedDF.final)[names(mergedDF.final) == 'SRR1767370'] <- 'PC3.3'

head(mergedDF.final)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# # Please denote which of the samples are interested in
# sampleofinterest <- "PC3"
# # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




# ##########################
# FILTER TO KEEP READS WITH NOT EXP UPC VALUES FOR PC3 AND/OR LNCaP
# ##########################
# should filter on UPC here and later output the actual normalised counts (CPM)

# first LNCaP and PC3
# second LNCaP only
# third PC3 alone
class(mergedDF.final$NP.1) # numeric

# UPC of 0.5 is the filter/expression threshold choice of the SCAN.UPC in PNAS



# PC3 or PC3-and-LNCaP have UPC <0.5
if(sampleofinterest == "PC3") { 
  cat("your sample is PC3")
  filterUPC <- which(
    (mergedDF.final$NP.1 >= 0.5 & mergedDF.final$NP.2 >= 0.5 & mergedDF.final$NP.3 >= 0.5 & mergedDF.final$LNCaP.1 < 0.5 & mergedDF.final$LNCaP.2 < 0.5 & mergedDF.final$LNCaP.3 < 0.5 & mergedDF.final$PC3.1 < 0.5 & mergedDF.final$PC3.1 < 0.5 & mergedDF.final$PC3.1 < 0.5)  
    
    | (mergedDF.final$NP.1 >= 0.5 & mergedDF.final$NP.2 >= 0.5 & mergedDF.final$NP.3 >= 0.5 & mergedDF.final$LNCaP.1 >= 0.5 & mergedDF.final$LNCaP.2 >= 0.5 & mergedDF.final$LNCaP.3 >= 0.5 & mergedDF.final$PC3.1 < 0.5 & mergedDF.final$PC3.1 < 0.5 & mergedDF.final$PC3.1 < 0.5)  
  )
}

# LNCaP or PC3-and-LNCaP have UPC <0.5
if(sampleofinterest == "LNCaP") { 
  cat("your sample is LNCaP")
  filterUPC <- which(
    (mergedDF.final$NP.1 >= 0.5 & mergedDF.final$NP.2 >= 0.5 & mergedDF.final$NP.3 >= 0.5 & mergedDF.final$LNCaP.1 < 0.5 & mergedDF.final$LNCaP.2 < 0.5 & mergedDF.final$LNCaP.3 < 0.5 & mergedDF.final$PC3.1 < 0.5 & mergedDF.final$PC3.1 < 0.5 & mergedDF.final$PC3.1 < 0.5)  
    
    | (mergedDF.final$NP.1 >= 0.5 & mergedDF.final$NP.2 >= 0.5 & mergedDF.final$NP.3 >= 0.5 & mergedDF.final$LNCaP.1 < 0.5 & mergedDF.final$LNCaP.2 < 0.5 & mergedDF.final$LNCaP.3 < 0.5 & mergedDF.final$PC3.1 >= 0.5 & mergedDF.final$PC3.1 >= 0.5 & mergedDF.final$PC3.1 >= 0.5)  
  )
}



mergedDF.final.UPCfiltered <- mergedDF.final[filterUPC,]
head(mergedDF.final.UPCfiltered)

(mergedDF.final.UPCfiltered$symbol)

# View(mergedDF.final.UPCfiltered)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










