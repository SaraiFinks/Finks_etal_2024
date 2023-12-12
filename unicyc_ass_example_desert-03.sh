#!/bin/bash

#SBATCH --job-name=curto_assembly_d03
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=1G	
#SBATCH --error=curto_assembly_d03_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com   

module load unicycler/0.4.8

WORKDIR=/path/to/longandshortreads/files
OUTDIR=/path/to/output/hybdrid/assembly

THREAD=32

cd $WORKDIR

unicycler -1 DESERT-3_S89_L001_R1.fq.gz -2 DESERT-3_S89_L001_R2.fq.gz -l D03_output.fastq.gz --contamination lambda --threads ${THREAD} \
	-o $OUTDIR/D03_third


