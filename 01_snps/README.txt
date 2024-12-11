Scripts used to call and filter SNPs using an original bcftools pipeline. The scrits were run on the institutional cluster with qsub system. Please note that most of the scripts utilize local paths and modules specific to the cluster. 

#SHELL SCRIPTS USED TO MAP READS TO THE REFERENCE GENOME
00_fastq2finalBAM.sh      - script used to trim and map sequening reads to the reference genome. To run the scripts individual id, fastq_r1, fast1_r2, reference genome and bowtie2 seq reference need to be provided.

#SHELL SCRIPTS USED TO CALL SNPS AND FILTER THEM 
01_SNP_calling_wraper.sh  - defines all required files, some parameters and runs in parallel SNP calling in windows along the genome
SNP_calling_chunks.sh     - calls SNPs in windows 
SNP_merge.sh              - merges all produced vcf files into a single one
FilterAndNormalize.sh     - fitlers and normizlizes SNPs 
runBEAGLE_sepPOP.sh       - script used to impute genotopes for some of the analyses 
annotateSNPs.sh           - annotate SNPs with snpEff and protein coding annotation

#SHELL SCRIPT USED TO RUN ADIMUXTURE
02_Admixture.sh           - filters vcf and run admixture

#FST ANALYSES
03_Fst                    - scripts to calculate Fst and ZFst outliers and GO enrtichment analyses 

#PSMC ANALYSES
04_psmc                   - scripts used to run psmc analyses 

