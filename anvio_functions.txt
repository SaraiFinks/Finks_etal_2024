#!/bin/bash

#SBATCH --job-name=
#SBATCH -A 	
#SBATCH -p standard		 	
#SBATCH --nodes=1        	
#SBATCH --ntasks=1		 	
#SBATCH --cpus-per-task=32   
#SBATCH --mem-per-cpu=4G	
#SBATCH --error=anvio_error.txt
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=email.com

module load anaconda/2020.07
source activate anvio-7

WORKDIR=/path/to/sequences

cd $WORKDIR

for f in *.fasta;
do

    basefileID=${f%.fasta}

    anvi-script-reformat-fasta "$f" -o $OUTDIR/${basefileID}-fixed.fasta -l 0 --seq-type NT -r $OUTDIR/${basefileID}-deflines.tsv --simplify-names

done

for f in *.fasta;
do

    basefileID=${f%.fasta}

    anvi-gen-contigs-database -f "$f" --skip-mindful-splitting -o $OUTDIR/${basefileID}-contigs.db 

done

for f in *.db;
do

 anvi-run-ncbi-cogs -c "$f" --sensitive --num-threads 32

done

for f in *.db;
do

anvi-run-pfams -c  "$f" --num-threads 32

done

for f in *.db;
do

anvi-run-hmms -c  "$f" --num-threads 32

done

anvi-compute-functional-enrichment -p curto_chr-PAN.db \
                                   -g curto_chr-GENOMES.db \
                                   --annotation-source COG20_FUNCTION \
                                   --category-variable habitat \
                                   -o COG20_FUNCTION_enrichment_habitat.txt \
                                   --functional-occurrence-table-output functions-occurrence-frequency_habitat.txt

anvi-compute-functional-enrichment -p curto_plas-PAN.db \
                                   -g curto_plas-GENOMES.db \
                                   --annotation-source COG20_FUNCTION \
                                   --category-variable clade \
                                   -o COG20_FUNCTION_enrichment_clade.txt \
                                   --functional-occurrence-table-output functions-occurrence-frequency_clade.txt

anvi-gen-genomes-storage -e external-genomes-chr.txt \
                         -o curto_chr-GENOMES.db

anvi-gen-genomes-storage -e external-genomes-p.txt \
                         -o curto_plas-GENOMES.db

anvi-compute-genome-similarity -e external-genomes-chr.txt \
                               -o ANI \
                               -p curto_chr-PAN.db 

anvi-compute-genome-similarity -e external-genomes-p.txt \
                               -o ANI \
                               -p curto_plas-PAN.db 

anvi-get-sequences-for-gene-clusters -p curto_chr-PAN.db \
                                     -g curto_chr-GENOMES.db \
                                     --min-num-genomes-gene-cluster-occurs 39 \
                                     --max-num-genes-from-each-genome 1 \
                                     --concatenate-gene-clusters \
                                     --output-file curto_chr-SCGs.fa

anvi-get-sequences-for-gene-clusters -p Curtobacterium_plasmids-PAN.db \
                                     -g curto_plasmid-GENOMES.db \
                                     --min-num-genomes-gene-cluster-occurs 18 \
                                     --max-num-genes-from-each-genome 1 \
                                     --concatenate-gene-clusters \
                                     --output-file curto_plasmid-SCGs_v2.fa

iqtree -s curto_chr-SCGs-cleaned.fa \
       -nt 8 \
       -m WAG \
       -bb 1000

iqtree -s curto_plasmid-SCGs.fa \
       -nt 8 \
       -m WAG \
       -bb 1000

trimal -in curto_chr-SCGs_v2.fa \
       -out curto_chr-SCGs-cleaned.fa \
       -gt 0.50

trimal -in curto_plasmid-SCGs.fa \
       -out curto_plasmid-SCGs-cleaned.fa \
       -gt 0.1

anvi-pan-genome -g curto_chr-GENOMES.db \
                --project-name "curto_chr" \
                --exclude-partial-gene-calls \
                --num-threads 32 \
                --minbit 0.5 \
                --mcl-inflation 2 \
                --sensitive \
                --min-occurrence 1 \
                --enforce-hierarchical-clustering

anvi-pan-genome -g curto_plas-GENOMES.db \
                --project-name "curto_plas" \
                --exclude-partial-gene-calls \
                --num-threads 32 \
                --minbit 0.5 \
                --mcl-inflation 2 \
                --sensitive \
                --min-occurrence 1 \
                --enforce-hierarchical-clustering

for f in *.db;
do

    basefileID=${f%.db}
 
    anvi-export-table --table gene_functions -o ${basefileID}-cog.txt ${basefileID}.db

done

for f in *.db;
do

    anvi-export-table --table hmm_hits -o hmm_hits_table.txt $f

done



