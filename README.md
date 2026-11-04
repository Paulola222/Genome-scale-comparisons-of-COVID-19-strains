# Genome-Scale Comparative Genomics of COVID-19 Strains

## Overview
This project performs comparative genomic analysis across SARS-CoV-2 strains 
sampled from multiple countries, alongside related coronaviruses (SARS, MERS), 
to identify conserved and variant genomic regions. The analysis uses sequence 
comparison algorithms to generate all-pairs similarity matrices, allowing 
regions of high conservation and high variability to be distinguished across 
strains sampled from geographically distinct populations.

## Data
- SARS-CoV-2 genome sequences (FASTA) from Australia, Brazil, China, India, 
  and the USA
- Related coronavirus strains: SARS (2003, 2017) and MERS (2012, 2014)
- Human BRCA2 and opsin gene sequences, included as non-viral controls to 
  validate that similarity signal is specific to related viral genomes 
  rather than a methodological artifact

## Method
Sequence comparison is performed using longest-common-substring identification 
and depth-first/postorder traversal of the resulting comparison structure to 
extract conserved and variant genomic regions across all strain pairs.

## Output
- `similarity_results.csv` — pairwise genomic similarity scores across all 
  strains compared
- `lcs_results.csv` — longest common substring results per strain pair
- `dfs_output.txt`, `postorder_output.txt` — traversal outputs used to 
  identify conserved regions

## Author
Paul Ayodeji Ola
