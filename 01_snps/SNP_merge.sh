#nfiles=$(wc -l ChunksSBS/splitchr.txt | awk '{print $1}')

zgrep "^#" ChunksSBS/tmp.1.SNPs.vcf.gz > headers.vcf

for i in {1..30737}; do zgrep -v "^#" ChunksSBS/tmp.${i}.SNPs.vcf.gz; done > tmp.all.variants.vcf

cat headers.vcf tmp.all.variants.vcf > all.population.SBS.variants.vcf 

rm -f headers.vcf
rm -f tmp.all.variants.vcf

#CHECK OUTPUT AND REMOVE ALL CHUNKS IN THE DIR Chunks SBS
