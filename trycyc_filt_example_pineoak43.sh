#!/bin/bash

#SBATCH --job-name=P43_filt	
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=1G	
#SBATCH --error=P43_filt_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com   

module load anaconda/2020.07
source activate filtlong


WORKDIR=/path/to/longreads

cd $WORKDIR

filtlong --min_length 1000 --keep_percent 95 P43_nanopore.fastq.gz > P43_reads.fastq

