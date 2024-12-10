#!/bin/bash

#+++++++++++++++++++++++++
#+++++++++++++++++++++++++
#+++++++++++++++++++++++++

FQPATH=/users/fk/mkonczal/projects/SBS/data/trimReads/batch4/    #NEED TO BE SPECIFIED: PATH TO FASTQFILES

IDNAME=$1     #INDIVIDUAL ID    $1x

FASTQ_1=$2    #FASTQ R1 FILE    $2
FASTQ_2=$3    #FASTQ R2 FILE    $3

REFSEQ=$4     #FASTA ASSEMBLY   $4      IT SHOULD BE INDEXED (SAMTOOLS); GATK DICTIONARY, BOWTIE2 INDEXES SHOULD BE PREPARED
BOWTIE2REF=$5 #BOWTIE2 SEQ REF  $5


#+++OPTIONAL OR TO CHANGE+++++

JMEM_OPT="-Xmx10G"
NTHREADS=1

#+++++++++++++++++++++++++++


echo Starting at `date`

#=======DIR FOR SAMPLE===============

mkdir $IDNAME
cd $IDNAME


#====================================
#====TRIM READS WITH TRIMMOMATIC=====
#====================================
echo Starting trimming reads at `date`

cp /software/fk/el6.3/Trimmomatic-0.32/adapters/TruSeq3-PE.fa .
java -jar /software/fk/el6.3/Trimmomatic-0.32/trimmomatic-0.32.jar PE $FQPATH/$FASTQ_1 $FQPATH/$FASTQ_2 \
${FASTQ_1}_trimmed.fastq_P.gz ${FASTQ_1}_trimmed.fastq_U.gz \
${FASTQ_2}_trimmed.fastq_P.gz ${FASTQ_2}_trimmed.fastq_U.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:30

echo Finished trimming reads at `date`

#====================================
#=========MAPPING====================  
#====================================  
echo Starting mapping at `date`

module load Bowtie2/2.2.3-goolf-1.4.10-no-OFED

R1=${FASTQ_1}_trimmed.fastq_P.gz
R2=${FASTQ_2}_trimmed.fastq_P.gz

bowtie2 -p $NTHREADS -x $BOWTIE2REF -1 $R1 -2 $R2 -S out.sam
samtools view -bS out.sam > out.bam
samtools sort -o ${IDNAME}.sorted.bam out.bam
samtools index ${IDNAME}.sorted.bam
samtools rmdup ${IDNAME}.sorted.bam $IDNAME.sorted.rmdup.bam 
samtools index ${IDNAME}.sorted.rmdup.bam

#rm -f out.sam
#rm -f out.bam               
#rm -f ${IDNAME}.sorted.bam     #OPTIONAL YOU MAY WANT TO KEEP ONE OF THEM FOR SOME REASONS
#rm -f ${IDNAME}.sorted.bam.bai

echo Finished mapping at `date`

#====================================
#====RQ ID AND REALIGMENT============
#====================================

#======RQ REPLACEMENT================ 
echo Starting RQ replacement at `date`

module load Java/1.8.0_74
export _JAVA_OPTIONS=$JMEM_OPT

java -jar /users/fk/mkonczal/software/picard-tools-2.1.1/picard.jar AddOrReplaceReadGroups I=$IDNAME.sorted.rmdup.bam O=$IDNAME.rg.bam RGID=$IDNAME RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=1 >> picard.log 2>&1

BAM1=$IDNAME.rg.bam

rm -f $IDNAME.sorted.rmdup.bam    ##OPTIONAL YOU MAY WANT TO KEEP ONE OF THEM FOR SOME REASONS
rm -f $IDNAME.sorted.rmdup.bam.bai

echo Finished RQ replacement at `date`

#======INDEX BAM====================
echo Sataring indexing bam with rq field  at `date`

samtools index $BAM1

echo Finished indexing bam with rq field  at `date`


#======INDELS REALIGMENT=============
echo Sataring indels realigment  at `date`

java -jar /users/fk/mkonczal/software/GenomeAnalysisTK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REFSEQ -I $BAM1 -o $IDNAME.intervals
java -jar /users/fk/mkonczal/software/GenomeAnalysisTK/GenomeAnalysisTK.jar -T IndelRealigner -R $REFSEQ -I $BAM1 -targetIntervals $IDNAME.intervals -o $IDNAME.final.bam

rm -f $BAM1
rm -f $BAM1.bai

#================================================
echo Finished indels religment at `date`
echo Done at `date`
