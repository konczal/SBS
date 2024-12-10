##first annotation file and reference genome used to create a "SBSv1c" database for snpEff

vcftools --vcf all.population.SBS.variants.FNNS.recalcAC.vcf --remove-indels --recode --recode-info-all --stdout | python VCF4annotation.py > all.FNNS.recAC.onlySNPs.BEAGLE4annotaton.vcf

VCF=all.FNNS.recAC.onlySNPs.BEAGLE4annotaton.vcf
java -jar /users/fk/mkonczal/software/snpEff/snpEff.jar ann -canon -ud 2000 SBSv1c ${VCF} > ${VCF}.annotated.vcf


