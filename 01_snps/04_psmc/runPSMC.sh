#!/usr/bin/bash

#PSMC for Ep sample 
samtools mpileup -uf SBS_final.scaffolds.fasta C_pyg_Ep.realined.bam | bcftools call -c | perl filterVCFonQualitymincov10.pl /dev/stdin 25 20 35 | vcfutils.pl vcf2fq > C_pyg_ep_consensus.fq
/users/fk/mkonczal/software/psmc/psmc -N50 -t5 -r5 -p "4+30*2+4+6+10" -o C_pyg_ep.psmc C_pyg_ep.psmcfa 
seq 100 | xargs -P 6 -i /users/fk/mkonczal/software/psmc/psmc -N50 -t5 -r5 -b -p "4+30*2+4+6+10" -o round-{}.psmc split.psmcfa

#PSMC for 2nd high cov. sample
samtools mpileup -uf ../SBS_final.scaffolds.fasta bird1.final.bam | bcftools call -c | perl ../filterVCFonQualitymincov10.pl /dev/stdin 25 20 38 | vcfutils.pl vcf2fq > C_pyg_bird1_consensus.fq


#PSMC for "C_pyg_embr" sample:
samtools mpileup -uf ../SBS_final.scaffolds.fasta embr.final.bam | bcftools call -c | perl ../filterVCFonQualitymincov10.pl /dev/stdin 25 20 32 | vcfutils.pl vcf2fq > C_pyg_embr.fq
echo "/users/fk/mkonczal/software/psmc/psmc -N50 -t5 -r5 -b -p "4+30*2+4+6+10" -o round-${SGE_TASK_ID}.psmc split_C_pyg_embr.psmcfa" > psmcM.sh
qsub -N array_psmc_test -t 1-100 -l virtual_free=15G,h_rt=12:00:00 -q long-sl65,fk-el6 -cwd psmcM.sh


#PSMC for C_ruf:
samtools mpileup -uf SBS_final.scaffolds.fasta C_ruf_02.realined.bam | bcftools call -c | perl filterVCFonQualitymincov10.pl /dev/stdin 25 20 26 | vcfutils.pl vcf2fq > C_ruf_consensus.fq
/users/fk/mkonczal/software/psmc/psmc -N50 -t5 -r5 -b -p "4+30*2+4+6+10" -o C_ruf.psmc C_ruf_consensus.psmcfa
echo "/users/fk/mkonczal/software/psmc/psmc -N50 -t5 -r5 -b -p "4+30*2+4+6+10" -o round-${SGE_TASK_ID}.psmc split.psmca" > psmcM_Cruf.sh
qsub -N array_psmc_test -t 1-100 -l virtual_free=15G,h_rt=12:00:00 -q long-sl65,fk-el6 -cwd psmcM_Cruf.sh 
