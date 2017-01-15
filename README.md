# PCaCellLineWGS_MS
Supplemental information, code, and data for the MS: 
Inge Seim, Penny L. Jeffery, Patrick B. Thomas, Colleen C. Nelson, Lisa K. Chopin. *Whole-genome sequence of the metastatic PC3 and LNCaP human prostate cancer cell lines*.

*Summary*: 
See individual scripts for additional information. Genome assemblies (in FASTA format available at ZENODO) and a BLAST server can be found at [ghrelinlab.org](http://ghrelinlab.org).  

### *de novo* genome assembly
`1A-sga-assembly.sh`: *de novo* genome assembly using [SGA](http://genome.cshlp.org/content/22/3/549.long). Formatted for usage on a HPC usage, but requires <150GB RAM and will work on modern desktops.

`1B-sga-gapfill.sh`: used to gapfill (add N's) a genome (in FASTA format, e.g. generated by `1A-sga-assembly.sh`)


### Variant detection
`2A-reads-vs-reference-genome.sh`: align FASTQ reads to a reference genome (FASTA) and call SNV and indels. Generates BAM and VCF files.

'2B-remove-common-variants.sh`: remove variants present that are likely common germline variants.

'2C-compare-and-annotate-sample-variants.sh`: compare variants present in two or more samples, filter using SnpSift and annodate using SnpEff.

3A-CNV-gene-loss-analysis.R
Requires indexed BAM files generated by `2A-reads-vs-reference-genome.sh` and the R package [cn.mops](http://nar.oxfordjournals.org/content/40/9/e69) (Copy Number estimation by a Mixture Of PoissonS) to identify copy number variation (CNV). Requires a BAM file generated from a 'normal' sample, here normal prostate. Calls `3B-load.RNAseq.R`, `biomaRt.look.up.R`, and `3C-merge.RNAseq.and.CNV.R` to generate a tab delimiated file of genes likely lost (CNV0, homozygous deletion events), including UPC values generated by the R package [SCAN.UPC](http://www.pnas.org/content/110/44/17778.long). SCAN.UPC outputs standardised expression values (UPC value), ranging from 0 to 1, which indicate whether a gene is actively transcribed in a sample of interest: higher values indicate that a gene is ‘active’. UPC scores are platform-independent and allow cross-experimental and cross-platform integration.

`4-cBioPortal.analysis.R`: Interrogates data downloaded from [cBioPortal](http://www.cbioportal.org/data_sets.jsp). Generates pie charts and performs survival analysis of genes stratified by copy number variation status (here CNV data on Capicua, *CIC*).