#!/bin/bash

#SBATCH --job-name=fastqc_curtos	
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=8   
#SBATCH --mem-per-cpu=1G	
#SBATCH --error=fastqc_error.txt 

module load fastqc/0.11.9

WORKDIR=/path/to/reads/
OUTDIR=/path/to/reads/

THREAD=8

cd $WORKDIR

for ffile in *.gz
do

	basefileID=${ffile%-READ1.fastq.gz}

	rfile=${basefileID}-READ2.fastq.gz

	fastqc $ffile $rfile --threads ${THREAD} \
	-o $OUTDIR
	
done