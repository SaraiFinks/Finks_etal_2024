#!/bin/bash

#SBATCH --job-name=P43_trycyc_polc
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=4G	
#SBATCH --error=P43_trycyc_polc_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com   

module load anaconda/2020.07
source activate masurca

WORKDIR=/path/to/trycycler/cluster_001

cd $WORKDIR

THREADS=32

polca.sh -a P43_polypolish.fasta -r "PINE-43_S65_L001_R1.fq.gz PINE-43_S65_L001_R2.fq.gz" -t ${THREADS} -m 4G
mv *.PolcaCorrected.fa P43_polypolish_polca.fasta