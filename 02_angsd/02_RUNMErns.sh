##this script is used to calculate effective number of syn and nonsyn sites as well as number of segregating sites
##as an innput we use angsd output generated across entire genome


##
ln -s ../04_nonsyn_syn_Sites/Scaffolds/
ln -s ../04_nonsyn_syn_Sites/E.pygmeus.evm.consensus.1c.gff3
ln -s ../04_nonsyn_syn_Sites/Scaffold_with_genes.txt
ln -s ../04_nonsyn_syn_Sites/GFF3/
ln -s ../04_nonsyn_syn_Sites/Codons/
ln -s ../04_nonsyn_syn_Sites/Proteins/
ln -s ../04_nonsyn_syn_Sites/Scaffold_names_dict.txt

mkdir VCF
mkdir EffectiveCodons

bcftools index ../03_angsd_SNPcalling/RNS_SNPs/RNS_AllSites.bcf

while read r1 r2 ; do
        bcftools view ../03_angsd_SNPcalling/RNS_SNPs/RNS_AllSites.bcf $r1 > VCF/$r2.vcf
        python ../Scripts/GetEffectiveCodons.py Codons/$r2.txt VCF/$r2.vcf > EffectiveCodons/${r2}_effCod.txt ; done < Scaffold_names_dict.txt

cd EffectiveCodons
find . -size 0 -delete
cd ..

##WARNING!! - it produces mutiple codon lines in case of multiple transcirpits... we need to correct for this (below we use only uniq codons)

mkdir EffectiveCodons_nonredundant

cd EffectiveCodons
for f in * ; do sort $f | uniq  > ../EffectiveCodons_nonredundant/${f}_nonred.txt ; done
cd ..

for f in EffectiveCodons_nonredundant/* ; do python ../Scripts/Count_NonsynSynSites.py $f ; done  > EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt

####Get VCF with SNPs within effective codons
mkdir VCF_SNPs

bcftools index ../03_angsd_SNPcalling/RNS_SNPs/RNS_SNPs.bcf

while read r1 r2; do
        bcftools view ../03_angsd_SNPcalling/RNS_SNPs/RNS_SNPs.bcf $r1 > VCF_SNPs/$r2.vcf ; done < Scaffold_names_dict.txt


ls EffectiveCodons/ | sed 's/_effCod.txt//g' > EffectiveCodons_Scaffolds.txt

mkdir VCF_SNPsINcodons
while read r1; do
  python ../Scripts/GetVCF.py EffectiveCodons_nonredundant/${r1}_effCod.txt_nonred.txt VCF_SNPs/${r1}.vcf > VCF_SNPsINcodons/${r1}_SNPs.vcf ; done < EffectiveCodons_Scaffolds.txt

grep "^#" VCF_SNPsINcodons/scaffold1_SNPs.vcf > header.vcf
grep -v "^#" VCF_SNPsINcodons/scaffold* > SNPs.vcf
cat header.vcf SNPs.vcf > SNPsInCodons_RNS.vcf
rm header.vcf SNPs.vcf

python ../Scripts/modifyVCF.py SNPsInCodons_RNS.vcf > SNPsInCodons_RNS_names.vcf
java -Xmx12g -jar ~/Software/SnpEff/snpEff/snpEff.jar SBS SNPsInCodons_RNS_names.vcf > SNPsInCodons_RNS_names.ANNOTATED.vcf

gzip -c SNPsInCodons_RNS_names.ANNOTATED.vcf > ../Results/SNPsInCodons_RNS_names.ANNOTATED.vcf.gz
gzip -c EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt > ../Results/EffectiveCodonsCounts_RNS.txt.gz
cd ../
##RNSS##
##THE FINAL FILES:
#SNPsInCodons_RNS_names.ANNOTATED.vcf                       - SNPs called in effective codons, you cna easly count them, also count number of nonsyn and syn, ans their are annotated
#EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt     - Effective number of codons per scaffold, you can easly sum them up and have all of them






