Scripts used to call and filter SNPs using an original bcftools pipeline. The scrits were run on the institutional cluster with qsub system.

00_fastq2finalBAM.sh      - script used to trim and map sequening reads to the reference genome 

01_SNP_calling_wraper.sh  - defines all required files, some parameters and runs in parallel SNP calling in windows along the genome
SNP_calling_chunks.sh     - calls SNPs in windows 
SNP_merge.sh              - merges all produced vcf files into a single one
FilterAndNormalize.sh     - fitlers and normizlizes SNPs 
runBEAGLE_sepPOP.sh       - script used to impute genotopes to some of the analyses 



