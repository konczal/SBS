# Population genomics of the critically endangered spoon-billed sandpiper

This repository contains code for population genomic analyses of the spoon-billed sandpiper and red-necked stint genomes.

## Data for the Analyses

- The required data, including the reference genome and its annotation, are available in the `00_Data` directory.
- An improved version of the reference genome is also available on [NCBI](https://www.ncbi.nlm.nih.gov/) (NCBI assembly: **ASM369795v1**), along with all sequencing reads used in the analyses (under the [NCBI Bioproject **PRJNA419629**](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA419629)).

## Provided Pipelines

Two analysis pipelines are included in this repository:

1. **BCFtools Pipeline**  
   Location: `01_snps` directory.

2. **ANGSD Pipeline**  
   Location: `02_angsd` directory.  
   Based on genotype likelihoods calculated with ANGSD.

Each pipeline directory contains a `README` file with detailed instructions and descriptions.

