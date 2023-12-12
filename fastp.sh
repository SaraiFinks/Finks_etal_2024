#!/bin/bash

#SBATCH --job-name=fastp 	
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=16   
#SBATCH --mem-per-cpu=1G	
#SBATCH --error=fastp_error.txt 


module load fastp/0.20.0

WORKDIR=/path/to/shortread/files
OUTDIR=/path/to/output/cleaned/reads

THREAD=16

cd $WORKDIR

for ffile in *.gz
do

	basefileID=${ffile%-READ1.fastq.gz}

	rfile=${basefileID}-READ2.fastq.gz

	fastp -i $ffile -I $rfile --average_qual 30 --length_required 20 --detect_adapter_for_pe -c \
	-o $OUTDIR/${basefileID}_R1.fq -O $OUTDIR/${basefileID}_R2.fq
	
done
