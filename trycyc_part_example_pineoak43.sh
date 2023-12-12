#!/bin/bash

#SBATCH --job-name=P43_trycyc_part
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=1G	
#SBATCH --error=P43_trycyc_part_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com   

module load anaconda/2020.07
source activate trycycler

WORKDIR=/path/to/partition/file

cd $WORKDIR

THREADS=32

trycycler partition --reads P43_reads.fastq --cluster_dirs trycycler/cluster_* --threads ${THREADS}
