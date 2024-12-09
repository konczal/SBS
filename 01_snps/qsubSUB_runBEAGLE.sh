##SC_NAME
##SCAFFOLD

WORK_DIR=/users/fk/mkonczal/projects/SBS/analyses/final_assembly/SNPcalling/SNPcalling_ALLDATA/BEAGLE_pop/
NTHREADS=5
CHRFILE=Scaffolds_SNPs.txt

####njobs=$(wc -l $CHRFILE  | awk '{print $1}')
##njobs=18403

njobs=18408
jobid=`qsub -t 1-${njobs} -pe smp $NTHREADS -e /dev/null -o /dev/null -cwd -l virtual_free=8G,h_rt=03:00:00 runBEAGLE_sepPOP.sh ${CHRFILE} ${NTHREADS}`
##jobid=`qsub -t 1-${njobs} -pe smp $NTHREADS -cwd -l virtual_free=7G,h_rt=03:00:00 runBEAGLE.sh ${CHRFILE} ${NTHREADS}`
##jobid=`qsub -t 1-3 -cwd -e /dev/null -o /dev/null -l virtual_free=15G,h_rt=03:00:00 runBEAGLE.sh ${CHRFILE} ${NTHREADS}`


