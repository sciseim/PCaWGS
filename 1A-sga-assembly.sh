# 1A-sga-assembly.sh
#!/bin/bash -l
#PBS -N sga-HPC
#PBS -l walltime=300:00:00
#PBS -l mem=150G
#PBS -l ncpus=12
#PBS -m bae
#PBS -j oe
module load bamtools/2.2.3 ;
module load bwa ;
module load samtools ;
module load abyss ;
module load bwa/0.7.7 ;

# set to you sga installation path
SGAPATH=/home/seim/bix/bin ; 
# The number of threads to use
CPU=12
# k-mer
Kmer=71 # same as in the Genome Research paper for the human genome
cd /home/seim/ProjectX/assembly/LNCaP/trimmed/PAIRED ;  # PAIRED-END FASTQ files
READ1=output_forward_paired.fq
READ2=output_reverse_paired.fq

# SGA error correction
# improved method from 'I want to assemble a very large genome. How do I efficiently index the reads?''
# https://github.com/jts/sga/wiki/FAQ

# #########################################
# 1. Preprocessing
#
# sga preprocess
# Prepare reads for assembly
# #########################################
$SGAPATH/sga preprocess --pe-mode 1 $READ1 $READ2 > mygenome.pp.fastq ;


# #########################################
# 2. Error Correction
# #########################################
# build index
# Build the FM-index for READS, which is a fasta or fastq file.
# This program is threaded (-t N).
# As the error corrector does not require the reverse BWT, suppress
# construction of the reversed index
# --no-reverse
# The general rule is that if you going to perform error correction, you should # run index with the "--no-reverse" flag.
# If you are going to assemble the reads into contigs (using sga filter, sga fm-merge, sga overlap) then you should *not* use the no-reverse flag.
# It doesn't matter whether the reads are paired or not when choosing the flags.
$SGAPATH/sga index -a ropebwt -t $CPU --no-reverse mygenome.pp.fastq ;

# #########################################
# Perform error correction on READS file. Overlap and kmer-based correction algorithms are implemented. By default, a k-mer based correction is performed.
# #########################################
$SGAPATH/sga correct -k $Kmer --learn -d 256 -t $CPU -p mygenome.pp --metrics=correct.metrics mygenome.pp.fastq > mygenome_correct.log ;


# #########################################
# Index the corrected data.
# #########################################
$SGAPATH/sga index -a ropebwt -t $CPU mygenome.pp.ec.fa  ;
$SGAPATH/sga filter -x 2 -t $CPU mygenome.pp.ec.fa  ;

# #########################################
# Merge together reads that can be unambiguously assembled
# #########################################
$SGAPATH/sga fm-merge -m 65 -t $CPU mygenome.pp.ec.filter.pass.fa

# #########################################
# Index the merged sequences
# #########################################
$SGAPATH/sga index -d 2000000 -t $CPU mygenome.pp.ec.filter.pass.merged.fa ;


# #########################################
# Use overlap to construct the string graph of the merged reads
# #########################################
$SGAPATH/sga overlap -t $CPU mygenome.pp.ec.filter.pass.merged.fa

# #########################################
# Perform the contig assembly
# #########################################
# The minimum overlap value used is 77. Aggressive variant removal parameters are chosen.
# Small repeats at the ends of reads are resolved using the -r 10 parameter. The minimum
# branch length for the trimming algorithm is 200 bp.
$SGAPATH/sga assemble -m 77 -d 0.4 -g 0.1 -r 10 -l 200 mygenome.pp.ec.filter.pass.merged.asqg.gz # WORKED


# The longest contig has the same length in all cases but N50 increases with the increase in the value of the -m parameter specified for the assemble command. Therefore, we will use the -m
# 75 and -m 77 assemblies.

# #########################################
# SCAFFOLDING
# #########################################
# Scaffolds consist of overlapping contigs separated by gaps of known length

# Change to a directory to hold the scaffolded assembly
cd ..
mkdir pe; 
cd pe;

# should be ../PAIRED
# Symlink the assembly from step 9
ln -s ../PAIRED/*-contigs.fa final-contigs.fa
ln -s ../PAIRED/*-graph.asqg.gz final-graph.asqg.gz

#
# Step 10: Index the contigs with bwa
bwa index -a bwtsw final-contigs.fa

# Symlink the raw reads files
ls ../PAIRED/*.fq | xargs -i ln -s {}

# skip sga-align
# bwa-mem here
# BWA-MEM is a new alignment algorithm for aligning sequence reads or long query sequences against a large reference genome such as human. It automatically chooses between local and end-to-end alignments, supports paired-end reads and performs chimeric alignment. The algorithm is robust to sequencing errors and applicable to a wide range of sequence lengths from 70bp to a few megabases. For mapping 100bp sequences, BWA-MEM shows better performance than several state-of-art read aligners to date.
# # https://github.com/jts/sga/issues/121
bwa mem -t $CPU ./final-contigs.fa ./$READ1 ./$READ2 | samtools view -F2304 -b -o reads.bam -