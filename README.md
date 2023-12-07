# SARS-CoV-2 GWAS Pipeline

![alt text](https://github.com/elder099/SARS-CoV-2-GWAS/blob/main/SC2_GWAS_Figure.png)


## **Step 1: GWAS_IGED_CR_JOIN_Filtering.R**
### Input:
    - Integrated Genomic Epi Dataset (IGED_Fill_in_R.csv)
    - IGED merged to phenotype data (IGED_CR_Merge_Mask.csv)
### Process:
    - Merge input datasets to get total # of identifiable entries w/phenotype
    - Filter data to only include IDable entries (w/EPI_ISL_###)
### Outputs:
    - Merged and filtered data (IGED_Full_Phenotype.csv)
___
## **Step 2: Match_and_Download.R**
### Input:
    - Output from Step 1, all IDable entries with 1/0 for Phen1 (IGED_Full_Phenotype.csv)
### Process:
    - Perform 1:3 case-control matching on a demographic using matchit
    - Extract gisaid entries in a loop using accessions from matched data
    - Deduplicate gisaid extracts by GISAID Acc #
    - (Insert % ref coverage filtering here)
### Outputs:
    - Demographic 1:3 case-control matched data (Matched_Dem2_Metadata.csv)
    - Sequences from gisaid extracts in fasta format (GISAID_Seqs.fasta)
___
## **Step 3: Pipelinev6.0.sh**
### Step 3a: MAFFT alignment
    - Input: Output from Step 2, gisaid sequences (GISAID_Seqs.fasta)
    - Outputs: Aligned fasta file (MAFFT_aligned.fasta)

### Step 3b: NextClade Tree Building
    - Input:
      - Database provided by NextClade (sars-cov-2_data)
      - Output from Step 3a, aligned fasta file (MAFFT_aligned.fasta)
    - Output:
      - Many outputs, but importantly Phylogenetic tree in json format (nextclade.auspice.json)
### Step 3c: Json-to-Newick Conversion
    - Input: Output from Step 3b, phylogenetic tree (nextclade.auspice.json)
    - Output: Phylogenetic tree now in newick format used by TreeWAS (Newicktree.nwk)
___
## **Step 4: TreeWASv6.0.R** (Pipelinev6.0.sh)
### Input:
    - Genetic dataset: Step 3a output, aligned fasta file (MAFFT_aligned.fasta)
    - Phenotype dataset: Step 1 output, matched demographic phenotype data (Matched_Dem2_Metadata)
    - Tree dataset: Step 3c output, phylogenetic tree (Newicktree.nwk)

### Process:
  - **Genetic dataset**:
      - Convert fasta file into binary snp matrix
      - Filter down to 0.5% minor allele frequency
      - Filter and format sample names
  - **Phenotype dataset**:
      - Format phenotype data to match genetic_dataset sample names
  - **Tree dataset**:
      - Format and filter tree tips to match genetic_dataset sample names
      - Filter tree tips to only include those in genetic and phenotype datasets
  - **TreeWAS:**
      - Run TreeWAS, correcting for lopsided phenotype distribution
      - Concatenate all significant SNPs for each score type, includes group counts

### Outputs:
    - Concatenated significant SNPs + group counts (SigSNPs.csv)
    - All SNPs, not all significant, no group counts (AllSNPs.csv)
___
## **Step 5: SNP Annotation**
### Inputs:
    - Nucleotide to AA converter, https://github.com/cednotsed/SARS-CoV-2-hookup (SARS-CoV-2_hookup_table_V3.csv)
    - Output of Step 4, list of SNPs in nucleotide format (AllSNPs.csv or SigSNPs.csv)
### Process:
    - Format SNP dataset to match vcf-like hookup table
    - Join SNP dataset and hookup table on position and nucleotide
    - Neatly format AA annotation
### Outputs:
    - All TreeWAS SNPs now in AA format (Annotated_SNPs.csv)
