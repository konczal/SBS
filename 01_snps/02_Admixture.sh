vcftools --gzvcf all.population.SBS.variants.FNNS.recalcAC.vcf.gz --keep SBS_RNS.ind.txt --remove-indels --minDP 6 --recode --recode-INFO-all --stdout | sed "s/|size[0-9]*//g" | sed "s|/nfs/users/fk/mkonczal/projects/SBS/analyses/final_assembly/Mapping/final_bams/||g" | sed "s|/users/fk/mkonczal/projects/SBS/analyses/final_assembly/Mapping/C_ruf/Sep2016/RNS_final_bams/||g" |sed "s|.realined.bam||g" | sed "s|.final.merged.bam||g" > all.population.filter.SBS.RNS.vcf

~/software/plink --vcf all.population.filter.SBS.RNS.vcf --allow-extra-chr --double-id --indep-pairwise 200 10 0.1 --maf 0.05
~/software/plink --vcf all.population.filter.SBS.RNS.vcf --extract plink.prune.in --make-bed --out prunedData --double-id --allow-extra-chr
~/software/plink --geno 0.2 --bfile prunedData --make-bed --out prunedData_filtered --double-id --allow-extra-chr
cp prunedData_filtered.bim prunedData_filtered.backup.bim
sed -i "s/scaffold//g" prunedData_filtered.bim
~/software/admixture_linux-1.3.0/admixture -j8 prunedData_filtered.bed 2
