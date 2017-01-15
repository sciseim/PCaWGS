# biomaRt look-up
ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

infotime=getBM(attributes = c("hgnc_symbol","description","entrezgene", "chromosome_name", "start_position", "end_position","gene_biotype"),filters = c("chromosomal_region"),values = list(chromosomal_region=filterlist), mart = ensembl)


# # works
# TESTfilterlist <- "Y:21820001:26395000"
# TESTinfotime=getBM(attributes = c("hgnc_symbol","description","entrezgene", "chromosome_name", "start_position", "end_position","gene_biotype"),filters = c("chromosomal_region"),values = list(chromosomal_region=TESTfilterlist), mart = ensembl)
# TESTinfotime
