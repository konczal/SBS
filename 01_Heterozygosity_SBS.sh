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


