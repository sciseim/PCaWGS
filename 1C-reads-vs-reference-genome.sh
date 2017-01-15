# 1C-reads-vs-reference-genome.sh

# The number of threads to use
CPU=10
# SAMPLE
SAMPLENAME=PC3
READ1=output_forward_paired.fq
READ2=output_reverse_paired.fq

REFERENCEGENOME=../GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

# Align the reads to the reference using BWA 
bwa mem -t$CPU $REFERENCEGENOME $READ1 $READ2 >$SAMPLENAME-reads.sam # 

# Call variants using bcftools.
# -@# is threads
samtools sort -@CPU -o $SAMPLENAME-reads.bam $SAMPLENAME-reads.sam
samtools index $SAMPLENAME-reads.bam

samtools mpileup -u -d 9999 -L 9999 -f $REFERENCEGENOME $SAMPLENAME-reads.bam >$SAMPLENAME-reads.bcf

bcftools call -c -v -Oz $SAMPLENAME-reads.bcf >$SAMPLENAME-reads.vcf.gz
bcftools index $SAMPLENAME-reads.vcf.gz
