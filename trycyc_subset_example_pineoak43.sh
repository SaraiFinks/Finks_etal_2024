#!/bin/bash

#SBATCH --job-name=trycyc_subset_P43	
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=1G	
#SBATCH --error=trycyc_P43_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com   

module load anaconda/2020.07
source activate trycycler

WORKDIR=/path/to/trycycler/subsample

cd $WORKDIR

trycycler subsample --reads P43_reads.fastq --out_dir read_subsets_P43
