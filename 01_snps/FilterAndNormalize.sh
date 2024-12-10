#filter low quality SNPs:
bcftools_filter -g 5 -O z -o all.population.SBS.variants.FILTERED.vcf.gz -i 'QUAL > 30 && DP > 40 && DP < 250 && MQSB > 0.0001' all.population.SBS.variants.vcf.gz

bcftools norm -d any -m + -o all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf all.population.SBS.variants.FILTERED.NORMALIZED.vcf.gz
bgzip all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf
tabix all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf
bcftools stats all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf.gz > all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf.gz.stats

mv ll.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf.gz ll.population.SBS.variants.FNNS.vcf.gz
