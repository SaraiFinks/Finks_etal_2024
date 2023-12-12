#!/bin/bash

#SBATCH --job-name=cazyid 	
#SBATCH -A 	
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=4G	
#SBATCH --error=canzyid_error.txt
#SBATCH --mail-user=email.com

module load anaconda/2022.05
source activate run_dbcan

WORKDIR=/path/to/files
OUTDIR=/path/to/output
DB=/path/to/dbCAN/database

cd $WORKDIR

for ffile in *.faa
do
	genomeID=$(echo $ffile)

	run_dbcan $ffile protein \
	--dbCANFile $DB/dbCAN.txt \
	--db_dir $DB \
	--tools all \
	--out_dir $OUTDIR/$genomeID

done
