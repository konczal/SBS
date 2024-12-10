#filter low quality SNPs:
bcftools_filter -g 5 -O z -o all.population.SBS.variants.FILTERED.vcf.gz -i 'QUAL > 30 && DP > 40 && DP < 250 && MQSB > 0.0001' all.population.SBS.variants.vcf.gz
#normalize:
bcftools_norm -f /users/fk/mkonczal/projects/SBS/assembly/final/SBS_final.scaffolds.fasta -D -o all.population.SBS.variants.FILTERED.NORMALIZEDrepNewSoft.vcf.gz -m + -O z all.population.SBS.variants.FILTERED.vcf.gz
#as I had some isses in the past with allele counts calculated by bcftools, I recalculate them here with my own python script:
zcat all.population.SBS.variants.FILTERED.NORMALIZEDrepNewSoft.vcf.gz | python ReCalculateAC.py | bgzip > all.population.SBS.variants.FNNS.recalcAC.vcf.gz
