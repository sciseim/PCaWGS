# 2C-compare-and-annotate-sample-variants.sh

# ############################################################### 
# COMPARE YOUR CELL LINE OF INTEREST
#
# # first merge into one
# Concatenate or combine VCF/BCF files. All source files must have the same sample columns appearing in the same order. Can # be used, for example, to concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel VCF into one. The # input files must be sorted by chr and position. The files must be given in the correct order to produce sorted VCF on output unless the -a, --allow-overlaps option is specified.

LNCAP1KVCF=/media/sciseim/mogwai/PCaWGS_MS/variants/LNCaP-reads/split_by_chr/1kg-filtered
PC31KVCF=/media/sciseim/mrwing/PC3-assembly/SnpEff/VCF-testing/PC3-reads/split_by_chr/1kg_germline_output

bcftools concat $LNCAP1KVCF/*.vcf  >LNCaP-1k-filtered.vcf
bcftools concat $PC31KVCF/*.vcf  >PC3-1k-filtered.vcf

# IDENTIFY COMMON SNVS AND INDELS
bcftools isec  -p output PC3-1k-filtered.vcf.gz LNCaP-1k-filtered.vcf.gz
# obtain
# output/0000.vcf	for records private to	PC3-1k-filtered.vcf.gz
# output/0001.vcf	for records private to	LNCaP-1k-filtered.vcf.gz
# output/0002.vcf	for records from PC3-1k-filtered.vcf.gz shared by both	PC3-1k-filtered.vcf.gz LNCaP-1k-filtered.vcf.gz
# ############################################################### 


# ###################################################
# SnpSift filter before annotating with SnpEff
# ###################################################
mkdir SnpSiftOut_no_annotyet ;

JAVAMEMORY=45g

java -jar -Xmx$JAVAMEMORY /bix/snpEff/SnpSift.jar filter "( QUAL >= 200 && DP >= 30 )" 0000.vcf >>./SnpSiftOut_no_annotyet/0000.vcf.snpsift.vcf && \
java -jar -Xmx$JAVAMEMORY /bix/snpEff/SnpSift.jar filter "( QUAL >= 200 && DP >= 30 )" 0001.vcf >>./SnpSiftOut_no_annotyet/0001.vcf.snpsift.vcf && \
java -jar -Xmx$JAVAMEMORY /bix/snpEff/SnpSift.jar filter "( QUAL >= 200 && DP >= 30 )" 0002.vcf >>./SnpSiftOut_no_annotyet/0002.vcf.snpsift.vcf && \

cd SnpSiftOut_no_annotyet ;

# annotate
mkdir ./stats && mkdir ./SnpEffOut && \
for VCF in *.vcf ; do java -jar -Xmx$JAVAMEMORY /bix/snpEff/snpEff.jar eff hg38 -v $VCF >./SnpEffOut/$VCF.snpeff -stats ./stats/$VCF.html; done
# ############################################################### 