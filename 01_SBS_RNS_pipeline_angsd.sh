##Spoon-billed sanpiper (SBS) genome project
##ANGSD-based analyses to calculate heterzygosity and density pf SNPs for Spoon-billed sandpiper (SBS) and for its sister species Red-necked stint (RNS)
## Mateusz Konczal 2021 (mateusz.konczal@amu.edu.pl)

##########################################
##ESTIMATE EXPECTED DEPTH PER INDIVIDUAL##
##########################################
mdkir 02_angsd_HET
cd 02_angsd_HET

##Required files:
GENOME=/media/raid/home/mkonczal/Projects/SBS_2021/00_data/genome/SBS_final.scaffolds.fasta  #reference genome
BAMfiles=All_bams.txt                                                                        #txt file with all bam files

##Output dir:
#creating dirs is commented in this script, mostly beasue it is not fully autmated script
mkdir  TestDPdistr 

#scaffold1-5 
~/Software/angsd/angsd/angsd -P 4 -bam $BAMfiles -ref $GENOME -out TestDPdistr/scaffold1.qc -r 'scaffold1|size10631746' \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500
~/Software/angsd/angsd/angsd -P 4 -bam $BAMfiles -ref $GENOME -out TestDPdistr/scaffold2.qc -r 'scaffold2|size10659794' \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500
~/Software/angsd/angsd/angsd -P 4 -bam $BAMfiles -ref $GENOME -out TestDPdistr/scaffold3.qc -r 'scaffold3|size10744440' \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500
~/Software/angsd/angsd/angsd -P 4 -bam $BAMfiles -ref $GENOME -out TestDPdistr/scaffold4.qc -r 'scaffold4|size9568156' \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500
~/Software/angsd/angsd/angsd -P 4 -bam $BAMfiles -ref $GENOME -out TestDPdistr/scaffold5.qc -r 'scaffold5|size11272743' \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500

##To get expected values I used this Scripts/getMaxMinPerSample.R
##It calculates mean Coverage per bam file, min and max thresholds (minCoverage = max(4, meanCoverage/2); maxCoverage = ifelse(MeanCoverage >20.0, MeanCoverage*2, MeanCoverage*3) )

##RESULT: Results/MaxDpPersample.txt; cols: ID, minCov, maxCov, path to bam file
#######################################################################################################################
#C_pyg_11        4       26      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_pyg_11.realined.bam
#C_pyg_13        4       23      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_pyg_13.realined.bam
#C_pyg_18        4       25      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_pyg_18.realined.bam
#C_pyg_22        6       35      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_pyg_22.realined.bam
#C_pyg_26        5       32      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_pyg_26.realined.bam
#C_pyg_28        5       31      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_pyg_28.realined.bam
#C_pyg_Ep        11      65      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_pyg_Ep.realined.bam
#C_ruf_01        4       13      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_ruf_01.final.merged.bam
#C_ruf_02        6       37      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_ruf_02.realined.bam
#C_ruf_03        4       8       /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_ruf_03.final.merged.bam
#C_ruf_05        4       15      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_ruf_05.final.merged.bam
#C_ruf_06        4       16      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_ruf_06.final.merged.bam
#C_ruf_07        4       12      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_ruf_07.final.merged.bam
#C_ruf_08        4       10      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_ruf_08.final.merged.bam
#C_ruf_09        4       14      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_ruf_09.final.merged.bam
#C_ruf_12        4       10      /media/raid/home/mkonczal/Projects/SBS_2021/00_data/bams/C_ruf_12.final.merged.bam
#######################################################################################################################

#The same analyses were made for "aditional files" - files that were excluded from some of the amalyses

###########################################
##GENOME WIDE HETERZYGOSITY IN WINDOWS#####
###########################################

MAXDP=../Results/MaxDpPersample.txt

while read r1 r2 r3 r4; do
        mkdir TestHetWholeGenome/$r1 
        BAMFILE=${r4}
        ~/Software/angsd/angsd/angsd -i $BAMFILE -ref $GENOME -anc $GENOME -dosaf 1 -gl 1 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20  -setMinDepth $r2 -setMaxDepth $r3  -doCounts 1 \
        -out TestHetWholeGenome/$r1/$r1
        ~/Software/angsd/angsd/misc/realSFS TestHetWholeGenome/$r1/$r1.saf.idx  >  TestHetWholeGenome/$r1/${r1}_est2.ml
        ~/Software/angsd/angsd/misc/realSFS saf2theta  TestHetWholeGenome/$r1/$r1.saf.idx -P 4 -fold 1 -outname TestHetWholeGenome/$r1/$r1 -sfs TestHetWholeGenome/$r1/${r1}_est2.ml
        ~/Software/angsd/angsd/misc/thetaStat do_stat  TestHetWholeGenome/$r1/$r1.thetas.idx  -win 20000 -step 20000  -outnames TestHetWholeGenome/$r1/$r1.ThetaWindows.gz ; done  < ${MAXDP}

#output with hets in windows, for example: 02_angsd_HET/TestHetWholeGenome/C_pyg_11/C_pyg_11.ThetaWindows.gz.pestPG
#the same analyses were made for other samples, not included in the main text

cd ..

###########################################
##SNP CALLING FOR SBS AND RNS##############
###########################################
mkdir 03_angsd_SNPcalling
cd 03_angsd_SNPcalling

##test coverage in the same way as above, but for joint bam files (within SBS and within RNS)
mkdir SBS_DPdistr RNS_DPdistr

#scaffold1 (done for scaffolds1..5, and separetly for SBS and RNS)
BAMS=Cpyg_bams.txt
~/Software/angsd/angsd/angsd -P 4 -bam $BAMS -ref $GENOME -out SBS_DPdistr/scaffold1.qc -r 'scaffold1|size10631746' \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500

###############################
##Coverage for SBS            #
#Mean=72.7 ; Min=36 ; Max=145 #
###                           #
##Coverage for RNS            #
#Mean=37.4 ; Min=19 ; Max=75  #
###############################

###SBS#####
mkdir SBS_SNPs

##Calling all sites (variable and invariable) meeting calling criteria
~/Software/angsd/angsd/angsd -P 4 -bam Cpyg_bams.txt -ref $GENOME -out SBS_SNPs/SBS_AllSites \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 4 -nThreads 50 \
        -doCounts 1 -dumpCounts 2 -setMinDepth 36 -setMaxDepth 145 \
        -GL 1  -doMaf 1 -doMajorMinor 1 -dopost 1 -doBcf 1

##Calling SNPs
~/Software/angsd/angsd/angsd -P 4 -bam Cpyg_bams.txt -ref $GENOME -out SBS_SNPs/SBS_SNPs \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 4 -nThreads 50 \
        -doCounts 1 -dumpCounts 2 -setMinDepth 36 -setMaxDepth 145 \
        -GL 1  -doMaf 1 -doMajorMinor 1 -dopost 1 -doBcf 1 -SNP_pval 1e-6

###RNS####
mkdir RNS_SNPs

##Calling all sites (variable and invariable) meeting calling criteria
~/Software/angsd/angsd/angsd -P 4 -bam Cruf_bams.txt -ref $GENOME -out RNS_SNPs/RNS_AllSites \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 4 -nThreads 50 \
        -doCounts 1 -dumpCounts 2 -setMinDepth 19 -setMaxDepth 75 \
        -GL 1  -doMaf 1 -doMajorMinor 1 -dopost 1 -doBcf 1

##Calling SNPs
~/Software/angsd/angsd/angsd -P 4 -bam Cruf_bams.txt -ref $GENOME -out RNS_SNPs/RNS_SNPs \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 4 -nThreads 50 \
        -doCounts 1 -dumpCounts 2 -setMinDepth 19 -setMaxDepth 75 \
        -GL 1  -doMaf 1 -doMajorMinor 1 -dopost 1 -doBcf 1 -SNP_pval 1e-6

cd ..

###########################################
##IDENTIFICATION NSYN-SYN SITES AND SNPS###
###########################################
mkdir 04_nonsyn_syn_Sites
cd 04_nonsyn_syn_Sites

#Creating required directories and files
mkdir Scaffolds Codons Proteins GFF3 VCF EffectiveCodons tmp

#split fasta into files
#./SBS_final.scaffolds_names.fasta file has modified seq. names that include only scaffold number without information about size
cd Scaffolds
python ../../Scripts/SplitFasta.py ../SBS_final.scaffolds_names.fasta
cd ..

##REQUIRED LINK TO THE ANNOTATION FILE IN GFF3 FORMAT
GFF=~/00_Data/Annotation/E.pygmeus.evm.consensus.1c.gff3
cut -f 1 $GFF | sort | uniq > Scaffold_with_genes.txt

#split gff into files
while read r1;  do awk -v SCAF="$r1" '$1 == SCAF && $3 == "CDS"' $GFF > GFF3/${r1}.gff3 ; done < Scaffold_with_genes.txt

#get codons per scaffold
while read r1 ; do
        python ../Scripts/gff_to_codonPositions.py GFF3/${r1}.gff3 Scaffolds/${r1}.fa > Codons/${r1}.txt
        mv GFF3/${r1}-proteins.fa Proteins/ ; done < Scaffold_with_genes.txt


#Get callable sites and limit codons only to those that are present in vcf file (all 3 bases)
#Scaffold_names_dict.txt - text file with names of scaffolds corresponding to two fasta files (that is redundant and should have been avoided), e.g.:
#scaffold1|size10631746  scaffold1

while read r1 r2 ; do
        bcftools view ../03_angsd_SNPcalling/SBS_SNPs/SBS_AllSites.bcf $r1 > VCF/$r2.vcf
        python ../Scripts/GetEffectiveCodons.py Codons/$r2.txt VCF/$r2.vcf > EffectiveCodons/${r2}_effCod.txt ; done < Scaffold_names_dict.txt


#remove empty files
cd EffectiveCodons
find . -size 0 -delete
cd ..

#Above commands produce mutiple codon lines in case of multiple transcirpits, so below it is reduced to nonredundant lines:
#mkdir EffectiveCodons_nonredundant
cd EffectiveCodons
for f in * ; do sort $f | uniq  > ../EffectiveCodons_nonredundant/${f}_nonred.txt ; done
cd ..

#Coun effective number of nonsyn and syn sites per scaffold
for f in EffectiveCodons_nonredundant/* ; do python ../Scripts/Count_NonsynSynSites.py $f ; done  > EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt
gzip -c EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt > ../Results/EffectiveCodonsCounts.txt.gz

###########04_nonsyn_syn_Sites####################################################
##EFFECTIVE CODONS FOR SBS: ./Results/EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt.gz
##################################################################################
####Get VCF with SNPs within effective codons
#mkdir VCF_SNPs VCF_SNPsINcodons

while read r1 r2; do
        bcftools view ../03_angsd_SNPcalling/SBS_SNPs/SBS_SNPs.bcf $r1 > VCF_SNPs/$r2.vcf ; done < Scaffold_names_dict.txt

while read r1; do
        python ../Scripts/GetVCF.py EffectiveCodons_nonredundant/${r1}_effCod.txt_nonred.txt VCF_SNPs/${r1}.vcf > VCF_SNPsINcodons/${r1}_SNPs.vcf ; done < EffectiveCodons_Scaffolds.txt

##Here we are getting all the SBS SNPs within effective codons:
grep "^#" VCF_SNPsINcodons/scaffold1_SNPs.vcf > header.vcf
grep -v "^#" VCF_SNPsINcodons/scaffold* > SNPs.vcf
cat header.vcf SNPs.vcf > SNPsInCodons_SBS.vcf
rm header.vcf SNPs.vcf

#Change scafofld names in the vcf file: 
python ../Scripts/modifyVCF.py SNPsInCodons_SBS.vcf > SNPsInCodons_SBS_names.vcf

##snpEff annotation:
java -Xmx12g -jar ~/Software/SnpEff/snpEff/snpEff.jar SBS SNPsInCodons_SBS_names.vcf > SNPsInCodons_SBS_names.ANNOTATED.vcf
gzip -c SNPsInCodons_SBS_names.ANNOTATED.vcf > ../Results/SNPsInCodons_SBS_names.ANNOTATED.vcf.gz

###########04_nonsyn_syn_Sites####################################################
##SNP RESULT FOR SBS: ./Results/SNPsInCodons_SBS_names.ANNOTATED.vcf.gz
##################################################################################

cd ..
###############################################
##IDENTIFICATION NSYN-SYN SITES AND SNPS RNS###
###############################################
mkdir 05_nonsyn_syn_Sites_RNS
cd 05_nonsyn_syn_Sites_RNS

#Steps are analagous to those given above for SBS
#All the analyses are summarized in the file 02_RUNMErns.sh
sh ../02_RUNMErns.sh

###########05_nonsyn_syn_Sites_RNS################################################
##RESULT FOR RNS: ./Results/EffectiveCodonsCounts.txt.gz
##RESULT FOR RNS: ./REsults/SNPsInCodons_RNS_names.ANNOTATED.vcf.gz
##################################################################################

cd ../
###########################################
##CALCULATING HETEROZGOSITY IN NSYN-SYN####
###########################################
mkdir 06_HETnsyn_SBS
cd 06_HETnsyn_SBS/

ln -s ../04_nonsyn_syn_Sites/ .
ls Codons > Scaffolds.txt

mkdir NonsynSites
while read r1; do
        python ../Scripts/GetEffNonsynPos.py Codons/$r1 > NonsynSites/$r1  ; done < Scaffolds.txt


ln -s ../02_angsd_HET/MaxDpPersample.txt .

##Heterozygosity is calculated only in fourfold degenerate sites (e.g. synonymous variation) and in non-degenerate sites (e.g. missense substitutions)
awk '$2 == 1.0 {print FILENAME "\t" $0}'  NonsynSites/* | sed 's|NonsynSites/||g' |  sed 's/.txt//g' | cut -f 1,2 > NonSynSites.Sites
awk '$2 == 0.0 {print FILENAME "\t" $0}'  NonsynSites/* | sed 's|NonsynSites/||g' |  sed 's/.txt//g' | cut -f 1,2 > SynSites.Sites

##Above files are redundant, so take uniq
sort NonSynSites.Sites | uniq | sort -k1,1 -k2,2n > NonSynSites.Uniq.Sites
sort SynSites.Sites | uniq | sort -k1,1 -k2,2n > SynSites.Uniq.Sites


#Modify names to longer forms and postiotions +1
python ../Scripts/modifyScafNames.py NonSynSites.Uniq.Sites > NonSynSites.Uniq.Names.Sites
python ../Scripts/modifyScafNames.py SynSites.Uniq.Sites > SynSites.Uniq.Names.Sites

#Index Sites
~/Software/angsd/angsd/angsd sites index NonSynSites.Uniq.Names.Sites
~/Software/angsd/angsd/angsd sites index SynSites.Uniq.Names.Sites


##RUN HET
mkdir NonSynHet SynHet

while read r1 r2 r3 r4; do
        mkdir NonSynHet/$r1
        BAMFILE=${r4}
        ~/Software/angsd/angsd/angsd -i $BAMFILE -ref $GENOME -anc $GENOME -dosaf 1 -gl 1 -Sites NonSynSites.Uniq.Names.Sites  \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20  -setMinDepth $r2 -setMaxDepth $r3  -doCounts 1 -doMajorMinor 1 \
        -out NonSynHet/$r1/$r1
        ~/Software/angsd/angsd/misc/realSFS  NonSynHet/$r1/$r1.saf.idx  >   NonSynHet/$r1/${r1}_est2.ml ; done  < MaxDpPersample.txt

while read r1 r2 r3 r4; do
        mkdir SynHet/$r1
        BAMFILE=${r4}
        ~/Software/angsd/angsd/angsd -i $BAMFILE -ref $GENOME -anc $GENOME -dosaf 1 -gl 1 -Sites SynSites.Uniq.Names.Sites  \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20  -setMinDepth $r2 -setMaxDepth $r3  -doCounts 1 -doMajorMinor 1 \
        -out SynHet/$r1/$r1
        ~/Software/angsd/angsd/misc/realSFS  SynHet/$r1/$r1.saf.idx  >   SynHet/$r1/${r1}_est2.ml ; done   < MaxDpPersample.txt

#_est2.ml files are SFS, which for a diploid sample gives number of heterzygous and both homozygous positions. It can be easly used to calculate heterozygosity

cd ..


##!!!!!!!!!!!!!!!!##
###FINAL FILES######
##!!!!!!!!!!!!!!!!##

####################
###NUMBER OF SNPs###
sh Scripts/GetStats_VCF.sh 04_nonsyn_syn_Sites/SNPsInCodons_SBS_names.ANNOTATED.vcf SBS
#SBS     misense_variants         18253
#SBS     synonymous_variants      39600
#SBS     stop_gained_variants     284
#SBS     start_lost_variants      53
#SBS     stop_lost_variants       33
#SBS     start_retained_variants  19
#SBS     stop_retained_variants   26
#SBS     other_variants   2

##COBMINED WITH STOPS AND STARTS
#SBS    nonsynonymous_variants  18623
#SBS    synonymous_varinats     39645

sh Scripts/GetStats_VCF.sh 05_nonsyn_syn_Sites_RNS/SNPsInCodons_RNS_names.ANNOTATED.vcf RNS
#RNS     misense_variants         37377
#RNS     synonymous_variants      94116
#RNS     stop_gained_variants     651
#RNS     start_lost_variants      157
#RNS     stop_lost_variants       75
#RNS     start_retained_variants  44
#RNS     stop_retained_variants   72
#RNS     other_variants   8

##COBMINED WITH STOPS AND STARTS
#RNS    nonsynonymous_variants  38304
#RNS    synonymous_varinats     94232

###########################################
##NUMBER OF MISSENSE AND SYNONYMOUS SITES##
awk  '{sum+=$1;} END{print sum;}' 04_nonsyn_syn_Sites/EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt
#SBS    missense_sites  1.57461e+07

awk  '{sum+=$2;} END{print sum;}' 04_nonsyn_syn_Sites/EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt
#SBS    synonymous_sites        4.57976e+06

awk '{sum+=$1;} END{print sum;}' 05_nonsyn_syn_Sites_RNS/EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt
#RNS    missense_sites  1.83896e+07

awk '{sum+=$2;} END{print sum;}' 05_nonsyn_syn_Sites_RNS/EffectiveCodons_nonredundant/EffectiveCodonsCounts.txt
#RNS    synonymous_sites        5.43179e+06


############################
###SNP DENSITIES############
##SBSnsyn=18623/1.57461e+07
##SBSsyn=39645/4.57976e+06
##RNSnsyn=38304/1.83896e+07
##RNSsyn=94232/5.43179e+06

#SBS    nonsyn_density  0.001183
#SBS    syn_desnity     0.008657
#RNS    nonsyn_density  0.002083
#RNS    syn_desnity     0.017348


#############################
##HETERZYGOSITY PER SAMPLE###

##Heterzygosity based on fourfold degenerate sites
while read r1 r2 ; do python Scripts/SummarizeHet.py 06_HETnsyn_SBS/SynHet/$r1/${r1}_est2.ml $r1 ; done < 06_HETnsyn_SBS/MaxDpPersample.txt
#Sample Heterozygosity  Number_heterozygous_sites
#C_pyg_11        0.002216        8446.07
#C_pyg_13        0.002374        8648.41
#C_pyg_18        0.002378        9145.94
#C_pyg_22        0.002444        9576.31
#C_pyg_26        0.002443        9570.5
#C_pyg_28        0.002134        8007.97
#C_pyg_Ep        0.002446        4093.49
#C_ruf_01        0.00473 12384.95
#C_ruf_02        0.004982        19561.64
#C_ruf_03        0.004476        7202.06
#C_ruf_05        0.004432        12221.42
#C_ruf_06        0.004774        15122.01
#C_ruf_07        0.00412 9034.59
#C_ruf_08        0.004325        7895.46
#C_ruf_09        0.004787        13660.61
#C_ruf_12        0.004157        7801.19

##Heterzygosity based on non-degenerate sites
while read r1 r2 ; do python Scripts/SummarizeHet.py 06_HETnsyn_SBS/NonSynHet/$r1/${r1}_est2.ml $r1 ; done < 06_HETnsyn_SBS/MaxDpPersample.txt
#C_pyg_11        0.000534        8418.81
#C_pyg_13        0.000572        8499.95
#C_pyg_18        0.000561        8752.34
#C_pyg_22        0.00055 8938.77
#C_pyg_26        0.000565        9030.38
#C_pyg_28        0.000477        7566.36
#C_pyg_Ep        0.000627        4390.05
#C_ruf_01        0.000978        10409.17
#C_ruf_02        0.000911        15039.55
#C_ruf_03        0.000817        5139.47
#C_ruf_05        0.000872        10070.43
#C_ruf_06        0.000857        11011.13
#C_ruf_07        0.000714        6377.04
#C_ruf_08        0.000745        5486.51
#C_ruf_09        0.000959        11138.13
#C_ruf_12        0.000735        5550.84



#####!!!!#############
#######################
###SOME EXTRA CHECKS###
##NUMBER OF STOP CODONS
for f in 04_nonsyn_syn_Sites/EffectiveCodons_nonredundant/scaffo* ; do awk '$5 == "*"' $f ; done | wc -l
## 13714
#
for f in 05_nonsyn_syn_Sites_RNS/EffectiveCodons_nonredundant/scaffo* ; do awk '$5 == "*"' $f ; done | wc -l
## 15553
#
##NUMBER OF START CODONS
for f in 04_nonsyn_syn_Sites/EffectiveCodons_nonredundant/scaffo* ; do awk '$5 == "M"' $f ; done | wc -l
## 149675
#
for f in 05_nonsyn_syn_Sites_RNS/EffectiveCodons_nonredundant/scaffo* ; do awk '$5 == "M"' $f ; done | wc -l
## 171724
#####################

grep "missense_" 04_nonsyn_syn_Sites/SNPsInCodons_SBS_names.ANNOTATED.vcf | cut -f 8 | cut -f 3 -d ";" | sed 's/AF=//g' >  04_nonsyn_syn_Sites/Freqs_SBS_missense.txt
grep -v "missense_" 04_nonsyn_syn_Sites/SNPsInCodons_SBS_names.ANNOTATED.vcf | grep "synonymous_" | cut -f 8 | cut -f 3 -d ";" | sed 's/AF=//g' >  04_nonsyn_syn_Sites/Freqs_SBS_synonymous.txt

grep "missense_" 05_nonsyn_syn_Sites_RNS/SNPsInCodons_RNS_names.ANNOTATED.vcf | cut -f 8 | cut -f 3 -d ";" | sed 's/AF=//g' >  05_nonsyn_syn_Sites_RNS/Freqs_RNS_missense.txt
grep -v "missense_" 05_nonsyn_syn_Sites_RNS/SNPsInCodons_RNS_names.ANNOTATED.vcf | grep "synonymous_" | cut -f 8 | cut -f 3 -d ";" | sed 's/AF=//g' >  05_nonsyn_syn_Sites_RNS/Freqs_RNS_synonymous.txt
##plot using 04_nonsyn_syn_Sites/plot_SFS.R

###Aditional SBS samples, plus Cruf and C minuta/subminuta estimates###

#!!!!#
##see 07_HETnsyn_SBSadditionalSamples/RUME.sh file##
#!!!!#


######################
##Heterzygosity and number of minor alleles in SBS and RNS
#####################

mkdir 11_RareSNPsPerGenome
cd 11_RareSNPsPerGenome
mkdir AncestralMajor  RNS  SBS
cd AncestralMajor

###Create reference with major alleles in fasta sequence
###It is done based on the previously generated SNPcalling that include information about allele frequeinces in each site

#Creat a link to such file for SBS
ln -s ../../03_angsd_SNPcalling/SBS_freqs/SBS_freqs.mafs.gz .
zcat SBS_freqs.mafs.gz > SBS_freqs.mafs

##use simple python script to replace sites
##this script has absolute specific paths to files, so each time it has to be updated, if used on other machines
python ReplaceNuclotFasta.py
rm SBS_freqs.mafs

##do the same for RNS
ln -s ../../03_angsd_SNPcalling/RNS_freqs/RNS_freqs.mafs.gz .
zcat RNS_freqs.mafs.gz > RNS_freqs.mafs
python ReplaceNuclotFasta_RNS.py
rm RNS_freqs.mafs


cd ../SBS
##C_pyg.txt and C_ruf.txt files have been created. Each include info about samples, path to bam files and ranges of coverage
##Links to NonsynSites and SynSites from ./../06_HETnsyn_SBS have been created. These files include infromation about Nonsyn and Syn sites, as described above

mkdir NonSynHet SynHet
mkdir ../RNS/NonSynHet ../RNS/SynHet

####Runing SFS estimation for each sample of SBS, where ancestral allele is major allele:
GENOME=/media/raid/home/mkonczal/Projects/SBS_2021/00_data/genome/SBS_final.scaffolds.fasta
ANC=/media/raid/home/mkonczal/Projects/SBS_2021/11_RareSNPsPerGenome/AncestralMajor/SBS_final.scaffolds_SBSmajor.fasta

while read r1 r2 r3 r4; do
        mkdir NonSynHet/$r1
        BAMFILE=${r4}
        ~/Software/angsd/angsd/angsd -i $BAMFILE -ref $GENOME -anc $ANC -dosaf 1 -gl 1 -Sites NonSynSites.Uniq.Names.Sites  \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20  -setMinDepth $r2 -setMaxDepth $r3  -doCounts 1 -doMajorMinor 1 \
        -out NonSynHet/$r1/$r1
        ~/Software/angsd/angsd/misc/realSFS  NonSynHet/$r1/$r1.saf.idx  >   NonSynHet/$r1/${r1}_est2.ml

        mkdir SynHet/$r1
        ~/Software/angsd/angsd/angsd -i $BAMFILE -ref $GENOME -anc $ANC -dosaf 1 -gl 1 -Sites SynSites.Uniq.Names.Sites \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20  -setMinDepth $r2 -setMaxDepth $r3  -doCounts 1 -doMajorMinor 1 \
        -out SynHet/$r1/$r1
        ~/Software/angsd/angsd/misc/realSFS  SynHet/$r1/$r1.saf.idx  >   SynHet/$r1/${r1}_est2.ml; done < C_pyg.txt


####Runing SFS estimation for each sample of RNS, where ancestral allele is major allele in RNS:
GENOME=/media/raid/home/mkonczal/Projects/SBS_2021/00_data/genome/SBS_final.scaffolds.fasta
ANC=/media/raid/home/mkonczal/Projects/SBS_2021/11_RareSNPsPerGenome/AncestralMajor/SBS_final.scaffolds_RNSmajor.fasta

while read r1 r2 r3 r4; do
        mkdir ../RNS/NonSynHet/$r1
        BAMFILE=${r4}
        ~/Software/angsd/angsd/angsd -i $BAMFILE -ref $GENOME -anc $ANC -dosaf 1 -gl 1 -Sites NonSynSites.Uniq.Names.Sites  \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20  -setMinDepth $r2 -setMaxDepth $r3  -doCounts 1 -doMajorMinor 1 \
        -out ../RNS/NonSynHet/$r1/$r1
        ~/Software/angsd/angsd/misc/realSFS  ../RNS/NonSynHet/$r1/$r1.saf.idx  >   ../RNS/NonSynHet/$r1/${r1}_est2.ml

        mkdir ../RNS/SynHet/$r1
        ~/Software/angsd/angsd/angsd -i $BAMFILE -ref $GENOME -anc $ANC -dosaf 1 -gl 1 -Sites SynSites.Uniq.Names.Sites \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20  -setMinDepth $r2 -setMaxDepth $r3  -doCounts 1 -doMajorMinor 1 \
        -out ../RNS/SynHet/$r1/$r1
        ~/Software/angsd/angsd/misc/realSFS  ../RNS/SynHet/$r1/$r1.saf.idx  >   ../RNS/SynHet/$r1/${r1}_est2.ml; done < C_ruf.txt


####SUMARIZE###
while read r1 r2 ; do python ../../Scripts/SummarizeHetAndDensity.py SynHet/$r1/${r1}_est2.ml $r1 $r2 synonynous ; done < C_pyg.txt
while read r1 r2 ; do python ../../Scripts/SummarizeHetAndDensity.py NonSynHet/$r1/${r1}_est2.ml $r1 $r2 "nonsyn" ; done < C_pyg.txt
while read r1 r2 ; do python ../../Scripts/SummarizeHetAndDensity.py ../RNS/SynHet/$r1/${r1}_est2.ml $r1 $r2 "synonynous" ; done < C_ruf.txt
while read r1 r2 ; do python ../../Scripts/SummarizeHetAndDensity.py ../RNS/NonSynHet/$r1/${r1}_est2.ml $r1 $r2 "nonsyn" ; done < C_ruf.txt

#Sample Min_cov Sites   Heterozygosity  Denisty #_Het_sites     #_Hom_rare_sites        #_all_sites
C_pyg_11        4       synonynous      0.00225 13995.51        8585.5  2705.0  3811764.0
C_pyg_13        4       synonynous      0.00241 14254.88        8781.16 2736.86 3643358.0
C_pyg_18        4       synonynous      0.00241 14911.94        9257.12 2827.41 3845485.0
C_pyg_22        6       synonynous      0.00247 14929.26        9661.77 2633.75 3918801.0
C_pyg_26        5       synonynous      0.00247 15394.64        9662.0  2866.32 3917345.0
C_pyg_28        5       synonynous      0.00222 10368.41        8317.39 1025.51 3753073.0
C_pyg_Ep        11      synonynous      0.00246 5758.18 4120.33 818.92  1673344.0
C_pyg_09        5       synonynous      0.00224 13155.27        8771.53 2191.87 3919355.0
C_pyg_12        4       synonynous      0.00233 14427.93        8455.48 2986.22 3635137.0
C_pyg_27        4       synonynous      0.00244 14676.06        8792.24 2941.91 3597163.0
C_pyg_29        5       synonynous      0.00234 13489.52        9457.87 2015.83 4049446.0
C_pyg_b1        19      synonynous      0.00255 13172.71        10934.91        1118.9  4294702.0
C_pyg_b2        17      synonynous      0.00279 14495.4 12259.48        1117.96 4391894.0
C_pyg_em        14      synonynous      0.00257 17099.02        11161.38        2968.82 4340219.0
C_pyg_11        4       nonsyn  0.00054 13366.73        8557.12 2404.81 15773408.0
C_pyg_13        4       nonsyn  0.00058 13536.74        8661.81 2437.46 14863606.0
C_pyg_18        4       nonsyn  0.00057 14046.29        8863.03 2591.63 15595615.0
C_pyg_22        6       nonsyn  0.00056 13787.62        9029.05 2379.29 16244868.0
C_pyg_26        5       nonsyn  0.00057 14149.78        9127.47 2511.15 15978692.0
C_pyg_28        5       nonsyn  0.0005  9660.44 7927.86 866.29  15855318.0
C_pyg_Ep        11      nonsyn  0.00063 6205.04 4419.44 892.8   7005743.0
C_pyg_09        5       nonsyn  0.00051 12095.35        8362.3  1866.52 16427866.0
C_pyg_12        4       nonsyn  0.00056 13344.93        8300.36 2522.29 14829928.0
C_pyg_27        4       nonsyn  0.0006  14061.0 8785.25 2637.87 14576396.0
C_pyg_29        5       nonsyn  0.00053 12336.84        8963.34 1686.75 16952071.0
C_pyg_b1        19      nonsyn  0.00057 11991.69        10046.54        972.57  17724355.0
C_pyg_b2        17      nonsyn  0.00064 13579.46        11620.48        979.49  18083679.0
C_pyg_em        14      nonsyn  0.00058 15916.7 10379.25        2768.72 17931701.0
C_ruf_01        4       synonynous      0.00485 19407.52        12691.99        3357.77 2618334.0
C_ruf_02        6       synonynous      0.00506 30337.01        19868.7 5234.16 3926212.0
C_ruf_03        4       synonynous      0.00464 12259.9 7465.78 2397.06 1608996.0
C_ruf_05        4       synonynous      0.00454 19732.88        12522.18        3605.35 2757709.0
C_ruf_06        4       synonynous      0.00492 24189.53        15578.97        4305.28 3167419.0
C_ruf_07        4       synonynous      0.00422 15123.56        9261.12 2931.22 2192844.0
C_ruf_08        4       synonynous      0.00441 13054.45        8050.92 2501.77 1825723.0
C_ruf_09        4       synonynous      0.0049  20960.37        13971.21        3494.58 2853749.0
C_ruf_12        4       synonynous      0.00427 12697.87        8021.5  2338.18 1876658.0
C_ruf_01        4       nonsyn  0.001   15696.79        10620.22        2538.28 10641098.0
C_ruf_02        6       nonsyn  0.00092 24248.22        15210.3 4518.96 16503287.0
C_ruf_03        4       nonsyn  0.00085 9094.78 5323.59 1885.6  6290911.0
C_ruf_05        4       nonsyn  0.0009  16135.93        10337.85        2899.04 11548033.0
C_ruf_06        4       nonsyn  0.00089 18749.88        11396.35        3676.77 12851589.0
C_ruf_07        4       nonsyn  0.00074 11378.14        6611.48 2383.33 8930464.0
C_ruf_08        4       nonsyn  0.00076 9431.13 5598.96 1916.09 7360713.0
C_ruf_09        4       nonsyn  0.00098 17012.99        11391.03        2810.98 11618947.0
C_ruf_12        4       nonsyn  0.00076 9485.91 5766.62 1859.64 7547763.0
