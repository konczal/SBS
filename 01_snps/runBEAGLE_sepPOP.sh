#!/bin/bash

###SGE_TASK_ID=1
module load Java/1.8.0_74
export _JAVA_OPTIONS=-Xmx3200M 

SC_FILE=$1     #<--File with scaffold names
NTHREADS=$2    #<--CPUs to use , must be at least 4 <- for beagle it is set by hard

##SCNAME="scaffold33"
##SCAFFOLD="scaffold33|size4241161"

SCNAME=`awk -v i=${SGE_TASK_ID} -v j=1 'FNR == i {print $j}' ${SC_FILE}`   #<--scaffold name without special character (to mkdir, create files etc)
SCAFFOLD=`awk -v i=${SGE_TASK_ID} -v j=2 'FNR == i {print $j}' ${SC_FILE}` #<--scaffold name without special character (to mkdir, create files etc)

WOR_DIR=/users/fk/mkonczal/projects/SBS/analyses/final_assembly/SNPcalling/SNPcalling_ALLDATA/BEAGLE_pop/
VCF=${WOR_DIR}/all.population.SBS.variants.FNNS.recalcAC.vcf.gz
IND_NAMES=${WOR_DIR}/all_samples_names.txt
GENOME=/nfs/users/fk/mkonczal/projects/SBS/assembly/final/SBS_final.scaffolds.fasta

SBSpop=${WOR_DIR}/SBS.txt
RNSpop=${WOR_DIR}/RNS.txt
RKNpop=${WOR_DIR}/RKN.txt
CMIpop=${WOR_DIR}/Cmin.txt
CSUpop=${WOR_DIR}/Csub.txt

mkdir ${WOR_DIR}/${SCNAME}
cd ${WOR_DIR}/${SCNAME}

tabix -h $VCF $SCAFFOLD | vcftools --vcf - --recode --recode-INFO-all --keep ${SBSpop} --stdout --remove-indels > ${SCNAME}.SBS.vcf
tabix -h $VCF $SCAFFOLD | vcftools --vcf - --recode --recode-INFO-all --keep ${RNSpop} --stdout --remove-indels > ${SCNAME}.RNS.vcf
tabix -h $VCF $SCAFFOLD | vcftools --vcf - --recode --recode-INFO-all --keep ${RKNpop} --stdout --remove-indels > ${SCNAME}.RKN.vcf
tabix -h $VCF $SCAFFOLD | vcftools --vcf - --recode --recode-INFO-all --keep ${CMIpop} --stdout --remove-indels > ${SCNAME}.CMI.vcf
tabix -h $VCF $SCAFFOLD | vcftools --vcf - --recode --recode-INFO-all --keep ${CSUpop} --stdout --remove-indels > ${SCNAME}.CSU.vcf

java -jar /users/fk/mkonczal/software/beagle/beagle.05Jul16.587.jar gl=${SCNAME}.SBS.vcf  niterations=8 nthreads=3 out=${SCNAME}.SBS.beagle.temp 
java -jar /users/fk/mkonczal/software/beagle/beagle.05Jul16.587.jar gl=${SCNAME}.RNS.vcf  niterations=8 nthreads=3 out=${SCNAME}.RNS.beagle.temp
java -jar /users/fk/mkonczal/software/beagle/beagle.05Jul16.587.jar gl=${SCNAME}.RKN.vcf  niterations=8 nthreads=3 out=${SCNAME}.RKN.beagle.temp
java -jar /users/fk/mkonczal/software/beagle/beagle.05Jul16.587.jar gl=${SCNAME}.CMI.vcf  niterations=8 nthreads=3 out=${SCNAME}.CMI.beagle.temp
java -jar /users/fk/mkonczal/software/beagle/beagle.05Jul16.587.jar gl=${SCNAME}.CSU.vcf  niterations=8 nthreads=3 out=${SCNAME}.CSU.beagle.temp

zcat ${SCNAME}.SBS.beagle.temp.vcf.gz | bgzip > ${SCNAME}.SBS.beagle.vcf.gz 
zcat ${SCNAME}.RNS.beagle.temp.vcf.gz | bgzip > ${SCNAME}.RNS.beagle.vcf.gz 
zcat ${SCNAME}.RKN.beagle.temp.vcf.gz | bgzip > ${SCNAME}.RKN.beagle.vcf.gz
zcat ${SCNAME}.CMI.beagle.temp.vcf.gz | bgzip > ${SCNAME}.CMI.beagle.vcf.gz 
zcat ${SCNAME}.CSU.beagle.temp.vcf.gz | bgzip > ${SCNAME}.CSU.beagle.vcf.gz

tabix ${SCNAME}.SBS.beagle.vcf.gz
tabix ${SCNAME}.RNS.beagle.vcf.gz
tabix ${SCNAME}.RKN.beagle.vcf.gz
tabix ${SCNAME}.CMI.beagle.vcf.gz
tabix ${SCNAME}.CSU.beagle.vcf.gz

bcftools merge -m all ${SCNAME}.SBS.beagle.vcf.gz ${SCNAME}.RNS.beagle.vcf.gz ${SCNAME}.CMI.beagle.vcf.gz ${SCNAME}.CSU.beagle.vcf.gz ${SCNAME}.RKN.beagle.vcf.gz | bgzip > ${SCNAME}.beagle.vcf.gz
tabix ${SCNAME}.beagle.vcf.gz


###GET SCAFFOLD_FASTA <<<<<<<<<<<<<<<<<<<<
##python ${WOR_DIR}/getFasta.py $GENOME  ${SCNAME}.fasta "${SCAFFOLD}"

##
## TO DO FILTER vcf files to have only SNPs!!! 
##

##while read r1 r2; do bcftools consensus -f ${SCNAME}.fasta -H 1 -s "${r2}" ${SCNAME}.beagle.vcf.gz > ${r1}.consensus.${SCNAME}.fasta; done < ${IND_NAMES}
##while read r1 r2; do sed -i "s/>.*/>${r1}_entire_${SCNAME}/g" ${r1}.consensus.${SCNAME}.fasta ; done < ${IND_NAMES}  ###<<<CHANGE NAMES IN HEADERS 
##cat *.consensus.${SCNAME}.fasta | gzip > ${SCNAME}.all.haplotypes.fasta.gz   
##while read r1 r2; do rm -f ${r1}.consensus.${SCNAME}.fasta ; done < ${IND_NAMES}


##vcftools --gzvcf ${SCNAME}.beagle.vcf.gz --recode --recode-INFO-all --remove-indels --out ${SCNAME}.beagle.noindels.vcf
##python ${WOR_DIR}/VCF2FASTA.py ${SCNAME}.beagle.noindels.vcf.recode.vcf ${IND_NAMES} ${SCNAME} | gzip > ${SCNAME}.beagle.snps2fasta.fasta.gz


############################################
##OUTPUTS TO MOVE TO OUT DIRs##
##${SCNAME}.beagle.vcf.gz                   #vcf
##${SCNAME}.beagle.vcf.gz.tbi               #vcf index
##${SCNAME}.all.haplotypes.fasta.gz         #fasta scaffolds with changed nucleotides for each ind (haplotype 1)
##${SCNAME}.beagle.snps2fasta.fasta.gz      #fasta merged only SNPs (hplotype 1) 
############################################

##mv ${SCNAME}.beagle.vcf.gz ../VCF/
##mv ${SCNAME}.beagle.vcf.gz.tbi ../VCF/
##mv ${SCNAME}.all.haplotypes.fasta.gz ../HAPLOTYPES/
##mv ${SCNAME}.beagle.snps2fasta.fasta.gz ../VCF2FASTA/

mv ${SCNAME}.beagle.vcf.gz* ../VCF_beaglePop/
cd ${WOR_DIR}
rm -fr ${SCNAME}
