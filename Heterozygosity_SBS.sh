##Spoon-billed sanpiper (SBS) genome project
##ANGSD-based analyses to calculate heterzygosity and density pf SNPs for SBS and for its sister species Red-necked stint (RNS)
##(c) Mateusz Konczal 2021 (mateusz.konczal@amu.edu.pl)


##########################################
##ESTIMATE EXPECTED DEPTH PER INDIVIDUAL##
##########################################
cd 02_angsd_HET

##Required files:
GENOME=/media/raid/home/mkonczal/Projects/SBS_2021/00_data/genome/SBS_final.scaffolds.fasta  #reference genome
BAMfiles=All_bams.txt                                                                        #txt file with all bam files

##Output dir:
mkdir  TestDPdistr #creating dirs is commented in this script, mostly beasue it is not fully autmated script

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


