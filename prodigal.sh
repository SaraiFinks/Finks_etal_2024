#!/bin/bash

#SBATCH --job-name=prodigal	
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=1G	
#SBATCH --error=prodigal_error.txt

module load prodigal/2.6.3

WORKDIR=/path/to/fasta/assemblies
OUTDIR=/path/to/output/protein/seqs

cd $WORKDIR

for ffile in *.fasta
do

	genomeID=$(echo $ffile)

	prodigal -i $ffile -o $OUTDIR/$genomeID.genes -a $OUTDIR/$genomeID.faa
	
done