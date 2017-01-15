# 1B-sga-gapfill.sh
#!/bin/bash -l
#PBS -N gapfill
#PBS -l walltime=12:00:00
#PBS -l mem=60G
#PBS -l ncpus=1
#PBS -m bae
#PBS -j oe
module load bamtools/2.2.3 ;
module load bwa ;
module load samtools ;
module load abyss ;
module load bwa/0.7.7 ;

# The number of threads to use
CPU=1

# go to directory with your initial FASTA files
cd /home/seim/ProjectX/assembly/LNCaP/trimmed/PAIRED ;

time sga gapfill -o ../pe/scaffolds.m200.build1_gapfilled.fa --prefix=mygenome.pp.ec.filter.pass ../pe/scaffolds.m200.build1.fa

