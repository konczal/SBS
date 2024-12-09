bcftools norm -d any -m + -o all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf all.population.SBS.variants.FILTERED.NORMALIZED.vcf.gz
bgzip all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf
tabix all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf
bcftools stats all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf.gz > all.population.SBS.variants.FILTERED.NORMALIZED_NORMALIZED2.vcf.gz.stats
