VCF=all.FNNS.recAC.onlySNPs.BEAGLE.vcf.gz

#between species
vcftools --gzvcf ${VCF} --out SBS_RNS_fst --weir-fst-pop RNS.txt --weir-fst-pop SBS.txt --fst-window-size 25000 --fst-window-size 25000

#between pops in SBS
vcftools --gzvcf ${VCF} --fst-window-size 25000 --weir-fst-pop BelyakaSplit.txt --weir-fst-pop Meinopylgino.txt --keep BelyakaSplit.txt --keep Meinopylgino.txt --out Belyaka_Meinopylgino_fst

#between pops in RNS
vcftools --gzvcf ${VCF} --fst-window-size 25000 --weir-fst-pop RuskayaKoshka.txt --weir-fst-pop MeinopylginoRNS.txt --keep RuskayaKoshka.txt --keep MeinopylginoRNS.txt --out RuskayaKoshka_Meinopylgino_RNS_fst



