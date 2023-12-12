#!/bin/bash

#SBATCH --job-name=P43_trycyc_poly
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=4G	
#SBATCH --error=P43_trycyc_poly_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com   

module load anaconda/2020.07
source activate polypolish

WORKDIR=/path/to/cluster_001

cd $WORKDIR

THREADS=32

bwa index consensus.fasta
bwa mem -t ${THREADS} -a 8_medaka.fasta P43_R1.fastq.gz > alignments_1.sam
bwa mem -t ${THREADS} -a 8_medaka.fasta P43_R2.fastq.gz > alignments_2.sam
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish 8_medaka.fasta filtered_1.sam filtered_2.sam > P43_polypolish.fasta
