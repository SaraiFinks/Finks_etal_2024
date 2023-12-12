#!/bin/bash

#SBATCH --job-name=P43_assembly	
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=4G	
#SBATCH --error=P43_assembly_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com  

module load flye/2.7.1
module load anaconda/2020.07
source activate raven

WORKDIR=/path/to/longreads

cd $WORKDIR

THREADS=32  # change as appropriate for your system
mkdir assemblies

flye --nano-raw sample_01.fastq -g 4m --threads "$THREADS" --out-dir assembly_01 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && rm -r assembly_01
flye --nano-raw sample_02.fastq -g 4m --threads "$THREADS" --out-dir assembly_02 && cp assembly_02/assembly.fasta assemblies/assembly_02.fasta && rm -r assembly_02
raven --threads "$THREADS" sample_03.fastq > assemblies/assembly_03.fasta && rm raven.cereal

flye --nano-raw sample_04.fastq -g 4m --threads "$THREADS" --out-dir assembly_04 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && rm -r assembly_04
flye --nano-raw sample_05.fastq -g 4m --threads "$THREADS" --out-dir assembly_05 && cp assembly_05/assembly.fasta assemblies/assembly_05.fasta && rm -r assembly_05
raven --threads "$THREADS" sample_06.fastq > assemblies/assembly_06.fasta && rm raven.cereal

flye --nano-raw sample_07.fastq -g 4m --threads "$THREADS" --out-dir assembly_07 && cp assembly_07/assembly.fasta assemblies/assembly_07.fasta && rm -r assembly_07
raven --threads "$THREADS" sample_08.fastq > assemblies/assembly_08.fasta && rm raven.cereal
raven --threads "$THREADS" sample_09.fastq > assemblies/assembly_09.fasta && rm raven.cereal

flye --nano-raw sample_10.fastq --threads "$THREADS" --out-dir assembly_10 && cp assembly_10/assembly.fasta assemblies/assembly_10.fasta && rm -r assembly_10
raven --threads "$THREADS" sample_11.fastq > assemblies/assembly_11.fasta && rm raven.cereal
raven --threads "$THREADS" sample_12.fastq > assemblies/assembly_12.fasta && rm raven.cereal
