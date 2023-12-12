#!/bin/bash

#SBATCH --job-name=read_coverage_shortreads
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=1G	
#SBATCH --error=coverage_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com   

module load bwa/0.7.8
module load samtools/1.10
module load anaconda/2020.07
source activate qualimap

INDEX=/path/to/assembly/fasta/file
READ1=/path/to/forward/reads
READ2=/path/to/reverse/reads
OUTDIR=/path/to/output/coverage/files

#create reference index for mapping
bwa index $INDEX

# paired-end mapping, general command structure, adjust to your case
bwa mem $INDEX $READ1/DESERT-3_S89_L001_R1.fq.gz $READ2/DESERT-3_S89_L001_R2.fq.gz > $OUTDIR/pD03a.sam 

cd $OUTDIR

#clean-up read pairing information and flags
samtools sort -n -O sam pD03a.sam | samtools fixmate -m -O bam - pD03a.fixmate.bam

#delete the old sam file it takes up a lot of space
rm pD03a.sam

#convert to bam file and sort
samtools sort -O bam -o pD03a.sorted.bam pD03a.fixmate.bam

rm pD03a.fixmate.bam

#generate mapping statistics - if considering SNP analysis follow https://genomics.sschmeier.com/ngs-mapping/ tutorial - remove duplicate reads.
samtools flagstat pD03a.sorted.bam

#get read depth at all positions in the reference genome (how many reads are overlapping the genomic position)
samtools depth pD03a.sorted.bam | gzip > pD03a.depth.txt.gz

#using QualiMap version 2.2.2 to generate a visual of coverage by read mapping
qualimap bamqc -bam pD03a.sorted.bam
