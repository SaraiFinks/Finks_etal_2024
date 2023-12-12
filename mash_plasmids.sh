#!/bin/bash

#SBATCH --job-name=mash_plasmids	
#SBATCH -A 	
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=4G	
#SBATCH --error=mash_plasmid_error.txt   

module load anaconda/2022.05
source activate mash

WORKDIR=/path/to/fasta/files

cd $WORKDIR

	mash sketch -o reference CP018784.1_p.fasta CP080396.1_p.fasta NZ_CP041260.1_p.fasta NZ_CP045288.2_p.fasta NZ_CP045289.2_p.fasta NZ_CP045290.2_p.fasta NZ_CP066342.1_p.fasta NZ_CP071884.1_p.fasta NZ_CP074440.1_p.fasta NZ_CP081962.1_p.fasta NZ_CP081963.1_p.fasta NZ_CP088077.1_p.fasta NZ_CP088078.1_p.fasta NZ_CP088079.1_p.fasta pD03a.fasta pD03b.fasta pD26a.fasta pD26b.fasta pD35.fasta pD40.fasta pG07.fasta pG54a.fasta pG54b.fasta pP20.fasta pS16.fasta pW21.fasta \
	
	mash info reference.msh \

	mash dist reference.msh CP018784.1_p.fasta CP080396.1_p.fasta NZ_CP041260.1_p.fasta NZ_CP045288.2_p.fasta NZ_CP045289.2_p.fasta NZ_CP045290.2_p.fasta NZ_CP066342.1_p.fasta NZ_CP071884.1_p.fasta NZ_CP074440.1_p.fasta NZ_CP081962.1_p.fasta NZ_CP081963.1_p.fasta NZ_CP088077.1_p.fasta NZ_CP088078.1_p.fasta NZ_CP088079.1_p.fasta pD03a.fasta pD03b.fasta pD26a.fasta pD26b.fasta pD35.fasta pD40.fasta pG07.fasta pG54a.fasta pG54b.fasta pP20.fasta pS16.fasta pW21.fasta > distances.tab



 

 
