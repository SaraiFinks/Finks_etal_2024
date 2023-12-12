#!/bin/bash

#SBATCH --job-name=P43_trycyc_msa
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=4G	
#SBATCH --error=P43_trycyc_msa_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com   

module load anaconda/2020.07
source activate trycycler

WORKDIR=/path/to/msa

cd $WORKDIR

THREADS=32

trycycler msa --cluster_dir trycycler/cluster_001 --threads ${THREADS}
