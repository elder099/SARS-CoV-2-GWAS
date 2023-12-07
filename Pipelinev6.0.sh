#!/bin/bash

#Conda
. ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate nextclade

GWASdir='~/CDPH/GWAS/GWASv6.0'
GENOMES=$GWASdir'/GISAID_Seqs.fasta'
ref_fasta=$GWASdir'/Reference.fasta'


echo 'Run MAFFT for alignment'
#Run MAFFT to align SARS-CoV-2 Genomes
mafft --auto --keeplength --addfragments $GENOMES --thread 10 $ref_fasta > $GWASdir/MAFFT_aligned.fasta

echo 'Run NextClade for tree building'
#Run Nextclade to generate phylogenetic tree
nextclade run \
   --input-dataset $GWASdir/sars-cov-2_data \
   --output-all=$GWASdir/output/ \
   $GWASdir/MAFFT_aligned.fasta
   #$GWASdir/MAFFT_aligned.fasta

echo 'Run json-to-newick to convert for TreeWAS'
python $GWASdir/json-to-newick.py

echo 'Run TreeWAS'
Rscript $GWASdir/TreeWASv6.0.R
