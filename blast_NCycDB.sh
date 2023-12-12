#!/bin/bash

#SBATCH --job-name=blast_ncycldb 	
#SBATCH -A 		
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=4G	
#SBATCH --error=blast_ncycldb_error.txt   

module load ncbi-blast/2.10.0

WORKDIR=/path/to/file
OUTDIR=/path/to/file/output
BLASTDB=/path/to/ncycl/database

cd $WORKDIR

THREAD=32


#input for this must be aa sequences
for f in *.faa
do

	genomename=${f%.faa}

	echo "processing $genomename"

	blastp -query $f -db $BLASTDB/ncycldb -evalue 1e-5 -num_threads $THREAD \
	-outfmt '6 qseqid sseqid pident length mismatch qstart qend sstart send evalue qcovs qcovhsp' \
	-out $OUTDIR/$genomename.txt -max_target_seqs 1

done


 

 
