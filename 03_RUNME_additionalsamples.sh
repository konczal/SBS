mkdir 07_HETnsyn_SBSadditionalSamples
cd 07_HETnsyn_SBSadditionalSamples

ln -s ../04_nonsyn_syn_Sites/Codons .
ls Codons > Scaffolds.txt

mkdir NonsynSites
while read r1; do
  python ../Scripts/GetEffNonsynPos.py Codons/$r1 > NonsynSites/$r1  ;done < Scaffolds.txt

ln -s ../06_HETnsyn_SBS/NonsynSites .
ln -s ../Results/MaxDpPersample_Additional.txt

ln -s ../06_HETnsyn_SBS/NonSynSites.Uniq.Sites* .
ln -s ../06_HETnsyn_SBS/SynSites.Uniq.Names.Sites* .


##RUN HET
mkdir NonSynHet
mkdir SynHet

GENOME=/media/raid/home/mkonczal/Projects/SBS_2021/00_data/genome/SBS_final.scaffolds.fasta

while read r1 r2 r3 r4; do
        mkdir NonSynHet/$r1
        BAMFILE=${r4}
        ~/Software/angsd/angsd/angsd -i $BAMFILE -ref $GENOME -anc $GENOME -dosaf 1 -gl 1 -Sites NonSynSites.Uniq.Names.Sites  \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20  -setMinDepth $r2 -setMaxDepth $r3  -doCounts 1 -doMajorMinor 1 \
        -out NonSynHet/$r1/$r1
        ~/Software/angsd/angsd/misc/realSFS  NonSynHet/$r1/$r1.saf.idx  >   NonSynHet/$r1/${r1}_est2.ml ; done  < MaxDpPersample_Additional.txt


        #~/Software/angsd/angsd/misc/realSFS saf2theta   NonSynHet/$r1/$r1.saf.idx -P 4 -fold 1 -outname  NonSynHet/$r1/$r1 -sfs  NonSynHet/$r1/${r1}_est2.ml
        #~/Software/angsd/angsd/misc/thetaStat do_stat   NonSynHet/$r1/$r1.thetas.idx  -win 1 -step 1  -outnames  NonSynHet/$r1/$r1.ThetaWindows.gz ; done  < MaxDpPersample.txt

while read r1 r2 r3 r4; do
        mkdir SynHet/$r1
        BAMFILE=${r4}
        ~/Software/angsd/angsd/angsd -i $BAMFILE -ref $GENOME -anc $GENOME -dosaf 1 -gl 1 -Sites SynSites.Uniq.Names.Sites  \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20  -setMinDepth $r2 -setMaxDepth $r3  -doCounts 1 -doMajorMinor 1 \
        -out SynHet/$r1/$r1
        ~/Software/angsd/angsd/misc/realSFS  SynHet/$r1/$r1.saf.idx  >   SynHet/$r1/${r1}_est2.ml ; done   < MaxDpPersample_Additional.txt


        #~/Software/angsd/angsd/misc/realSFS saf2theta   NonSynHet/$r1/$r1.saf.idx -P 4 -fold 1 -outname  NonSynHet/$r1/$r1 -sfs  NonSynHet/$r1/${r1}_est2.ml
        #~/Software/angsd/angsd/misc/thetaStat do_stat   NonSynHet/$r1/$r1.thetas.idx  -win 1 -step 1  -outnames  NonSynHet/$r1/$r1.ThetaWindows.gz ; done  < MaxDpPersample.txt


##Heterzygosity based on fourfold degenerate sites
while read r1 r2 ; do python ../Scripts/SummarizeHet.py SynHet/$r1/${r1}_est2.ml $r1 ; done < MaxDpPersample_Additional.txt

#Sample Heterozygosity  Number_heterozygous_sites
#C_can_04        0.004135        15024.62
#C_min_18        0.00508 2512.5
#C_pyg_09        0.002185        8563.72
#C_pyg_12        0.002298        8354.33
#C_pyg_27        0.00241 8669.82
#C_pyg_29        0.002286        9257.86
#C_pyg_b1        0.002516        10805.37
#C_pyg_b2        0.002787        12240.33
#C_pyg_em        0.002561        11113.87
#C_sub_11        0.003899        1555.92

##Heterzygosity based on non-degenerate sites
while read r1 r2 ; do python ../Scripts/SummarizeHet.py NonSynHet/$r1/${r1}_est2.ml $r1 ; done < MaxDpPersample_Additional.txt
#C_can_04        0.000635        9858.92
#C_min_18        0.000742        1454.81
#C_pyg_09        0.000496        8141.3
#C_pyg_12        0.000552        8188.69
#C_pyg_27        0.000592        8623.21
#C_pyg_29        0.000513        8698.49
#C_pyg_b1        0.000559        9901.04
#C_pyg_b2        0.000641        11585.69
#C_pyg_em        0.000576        10323.4
#C_sub_11        0.000537        866.07

