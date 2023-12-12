#!/bin/bash

#SBATCH --job-name=mash_chromosomes	
#SBATCH -A 	
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=4G	
#SBATCH --error=mash_chromosomes_error.txt   

module load anaconda/2022.05
source activate mash

WORKDIR=/path/to/fasta/files

cd $WORKDIR

	mash sketch -o reference CP018783.1_c.fasta \
	CP045287.2_c.fasta \
	CP080395.1_c.fasta \
	CP093378.1_c.fasta \
	D03.fasta  \
	D26.fasta \
	D35.fasta \
	D40.fasta \
	G07.fasta \
	G31.fasta \
	G32.fasta \
	G36.fasta \
	G54.fasta \
	G58.fasta \
	MCBA15004.fasta \
	MCBA15012.fasta \
	MMLR14002_corrected.fasta \
	NZ_CP009755.1_c.fasta \
	NZ_CP017580.1_c.fasta \
	NZ_CP027869.1_c.fasta \
	NZ_CP041259.1_c.fasta \
	NZ_CP066341.1_c.fasta \
	NZ_CP068987.1_c.fasta \
	NZ_CP071883.1_c.fasta \
	NZ_CP074439.1_c.fasta \
	NZ_CP076544.1_c.fasta \
	NZ_CP081964.1_c.fasta \
	NZ_CP083910.1_c.fasta \
	NZ_CP088076.1_c.fasta \
	P20.fasta \
	P43.fasta \
	S05.fasta \
	S07.fasta \
	S15.fasta \
	S16.fasta \
	W02.fasta \
	W21.fasta \
	W50.fasta \
	W52.fasta \
	
	mash info reference.msh \

	mash dist reference.msh  CP018783.1_c.fasta \
	CP045287.2_c.fasta \
	CP080395.1_c.fasta \
	CP093378.1_c.fasta \
	D03.fasta  \
	D26.fasta \
	D35.fasta \
	D40.fasta \
	G07.fasta \
	G31.fasta \
	G32.fasta \
	G36.fasta \
	G54.fasta \
	G58.fasta \
	MCBA15004.fasta \
	MCBA15012.fasta \
	MMLR14002_corrected.fasta \
	NZ_CP009755.1_c.fasta \
	NZ_CP017580.1_c.fasta \
	NZ_CP027869.1_c.fasta \
	NZ_CP041259.1_c.fasta \
	NZ_CP066341.1_c.fasta \
	NZ_CP068987.1_c.fasta \
	NZ_CP071883.1_c.fasta \
	NZ_CP074439.1_c.fasta \
	NZ_CP076544.1_c.fasta \
	NZ_CP081964.1_c.fasta \
	NZ_CP083910.1_c.fasta \
	NZ_CP088076.1_c.fasta \
	P20.fasta \
	P43.fasta \
	S05.fasta \
	S07.fasta \
	S15.fasta \
	S16.fasta \
	W02.fasta \
	W21.fasta \
	W50.fasta \
	W52.fasta > distances.tab


mash triangle reference.msh  CP018783.1_c.fasta \
	CP045287.2_c.fasta \
	CP080395.1_c.fasta \
	CP093378.1_c.fasta \
	D03.fasta  \
	D26.fasta \
	D35.fasta \
	D40.fasta \
	G07.fasta \
	G31.fasta \
	G32.fasta \
	G36.fasta \
	G54.fasta \
	G58.fasta \
	MCBA15004.fasta \
	MCBA15012.fasta \
	MMLR14002_corrected.fasta \
	NZ_CP009755.1_c.fasta \
	NZ_CP017580.1_c.fasta \
	NZ_CP027869.1_c.fasta \
	NZ_CP041259.1_c.fasta \
	NZ_CP066341.1_c.fasta \
	NZ_CP068987.1_c.fasta \
	NZ_CP071883.1_c.fasta \
	NZ_CP074439.1_c.fasta \
	NZ_CP076544.1_c.fasta \
	NZ_CP081964.1_c.fasta \
	NZ_CP083910.1_c.fasta \
	NZ_CP088076.1_c.fasta \
	P20.fasta \
	P43.fasta \
	S05.fasta \
	S07.fasta \
	S15.fasta \
	S16.fasta \
	W02.fasta \
	W21.fasta \
	W50.fasta \
	W52.fasta > distances2.tab


 

 
