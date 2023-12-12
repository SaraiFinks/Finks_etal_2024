#!/bin/bash

#SBATCH --job-name=read_coverage	
#SBATCH -A 	
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=8   
#SBATCH --mem-per-cpu=1G	
#SBATCH --error=coverage_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com   

module load minimap2/2.23
module load samtools/1.10
module load anaconda/2020.07
source activate qualimap

INDEX=/path/to/fasta/file/shortreads
READ=/path/to/fasta/file/longreads
OUTDIR=/path/to/output/coverage/files

cd $OUTDIR

#create reference index for mapping (makes a minimizer index for the reference before mapping)
minimap2 -d D03l.mmi $INDEX

# paired-end mapping, general command structure, adjust to your case
minimap2 -ax map-ont D03l.mmi $READ/D03_long.fastq.gz > $OUTDIR/D03l.sam 

#clean-up read pairing information and flags
samtools sort -n -O sam D03l.sam | samtools fixmate -m -O bam - D03l.fixmate.bam

#delete the old sam file it takes up a lot of space
rm D03l.sam

#convert to bam file and sort
samtools sort -O bam -o D03l.sorted.bam D03l.fixmate.bam

rm D03l.fixmate.bam

#generate mapping statistics - if considering SNP analysis follow https://genomics.sschmeier.com/ngs-mapping/ tutorial - remove duplicate reads.
samtools flagstat D03l.sorted.bam

#get read depth at all positions in the reference genome (how many reads are overlapping the genomic position)
samtools depth D03l.sorted.bam | gzip > D03l.depth.txt.gz

#using QualiMap version 2.2.2 to generate a visual of coverage by read mapping
qualimap bamqc -bam D03l.sorted.bam

