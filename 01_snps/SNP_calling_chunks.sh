#!/bin/bash
#$ -cwd
##$ -M mateusz.konczal@crg.es
##$ -m a
#$ -V




echo Starting at `date`

export PATH=/software/fk/el6.3/samtools-1.3/bin/:/software/fk/el6.3/bcftools-1.3/bin/:/software/fk/el6.3/htslib-1.3/bin/:$PATH

bam_list=$1
genome_ref=$2
splitfile=$3

chrsplit=$(awk 'NR=='$SGE_TASK_ID'{print $1}' $splitfile)

samtools mpileup -C50 -R -t DP,ADF,ADR -vuf $genome_ref -r "${chrsplit}" -b $bam_list -o tmp.${SGE_TASK_ID}.var.raw.vcf
bgzip tmp.${SGE_TASK_ID}.var.raw.vcf && tabix tmp.${SGE_TASK_ID}.var.raw.vcf.gz

bcftools call -f GQ -vmO v -o tmp.${SGE_TASK_ID}.SNPs.vcf tmp.${SGE_TASK_ID}.var.raw.vcf.gz
bgzip tmp.${SGE_TASK_ID}.SNPs.vcf && tabix tmp.${SGE_TASK_ID}.SNPs.vcf.gz

rm -f tmp.${SGE_TASK_ID}.var.raw.vcf.gz
rm -f tmp.${SGE_TASK_ID}.var.raw.vcf.gz.tbi

echo Finished $step at `date`
