#!/usr/bin/env python
"""Convert GlimmerHMM GFF3 gene predictions into codon sequences and their coordinates
This works with the GlimmerHMM GFF3 output format:
##gff-version 3
##sequence-region Contig5.15 1 47390
Contig5.15      GlimmerHMM      mRNA    323     325     .       +       .       ID=Contig5.15.path1.gene1;Name=Contig5.15.path1.gene1
Contig5.15      GlimmerHMM      CDS     323     325     .       +       0       ID=Contig5.15.cds1.1;Parent=Contig5.15.path1.gene1;Name=Contig5.15.path1.gene1;Note=final-exon
http://www.cbcb.umd.edu/software/GlimmerHMM/
Usage:
    glimmergff_to_proteins.py <glimmer gff3> <ref fasta>

The script has been improved by mkonczal, so now, it also prints out codons (positions, codon (if reverse strand it is rc to what is in the reference sequence) and translated aa
in such a way it can be used to count effective number of codons and effective number of nonsynonyous and synonyous sites

STDOUT should look be in the following format:
##Epyg1c007437T1
4576773 4576774 4576775 CAG     Q
4576776 4576777 4576778 CTA     L
4576779 4576780 4576781 CTG     L
4576782 4576783 4576784 ACC     T

Addtionally a file with aa fasta sequences should be created
"""

from __future__ import with_statement
import sys
import os
import operator
from functools import reduce

from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

from BCBio import GFF

def main(glimmer_file, ref_file):
    with open(ref_file) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

    base, ext = os.path.splitext(glimmer_file)
    out_file = "%s-proteins.fa" % base
    with open(out_file, "w") as out_handle:
        SeqIO.write(protein_recs(glimmer_file, ref_recs), out_handle, "fasta")

def protein_recs(glimmer_file, ref_recs):
    """Generate protein records from GlimmerHMM gene predictions.
    """
    with open(glimmer_file) as in_handle:
        for rec in glimmer_predictions(in_handle, ref_recs):
            for feature in rec.features:
                print("##" + feature.id)
                seq_exons = []
                counter = 0
                frame = 0

##Getcodons
                cds_locs = []
##
                lastPhase = 0


                if feature.sub_features == []:
                    #print("FUCK")
                    feature.sub_features = [feature]

                for cds in feature.sub_features:
                    #print(cds.qualifiers['phase'][0])
                    lastPhase =  int(cds.qualifiers['phase'][0])
                    if cds.strand == 1:
                        if counter == 0:
                            frame =  int(cds.qualifiers['phase'][0]) # that is only for the first exon; the rest stay the same
                            counter += 1
                    #print(list(cds))

##Getcodons
                        cds_locs.extend(list(cds))
                    if cds.strand == -1:
                        cds_locs[0:0] = list(cds)
##

                    seq_exons.append(rec.seq[
                        cds.location.nofuzzy_start:
                        cds.location.nofuzzy_end])
                gene_seq = Seq.Seq(str(reduce(operator.add, seq_exons, "")))
                #print(cds_locs)
                #gene_seq = gene_seq2[frame:]  #modified here to check impact of strand -> there might be shift between CDS and ORF according to strand in the first exon (that is easy on + strand, but does not work on the strand "-")
                if feature.strand == -1:
                    gene_seq2 = gene_seq.reverse_complement()
                    gene_seq = gene_seq2[lastPhase:]
                    frame = 0
                protein_seq = gene_seq[frame:].translate()  #modified here for test

##Getcodons
                if feature.strand == -1:
                    #print(cds_locs)
                    #cds_locs.reverse()
                    #print(cds_locs)
                    #break
                    cds_locsUp = cds_locs[lastPhase:]
                if feature.strand == 1:
                    cds_locsUp = cds_locs[frame:]

                Codons = [cds_locsUp[i:i + 3] for i in range(len(cds_locsUp))[::3]]
                #print(Codons)
                for codon in Codons:
                    if len(codon) < 3:
                        continue
                    if feature.strand == 1:
                        fCod = rec.seq[codon[0]] + rec.seq[codon[1]] + rec.seq[codon[2]]
                        fCodSeq = Seq.Seq(fCod)
                    if feature.strand == -1:
                        fCod = rec.seq[codon[0]] + rec.seq[codon[1]] + rec.seq[codon[2]]
                        fCodSeq = Seq.Seq(fCod).complement()
                    fCodT = fCodSeq.translate()
                    print(" ".join([str(x) for x in codon]) + "\t" + str(fCodSeq) + "\t" + str(fCodT))
                #print("###")
##

                yield SeqRecord(protein_seq, feature.qualifiers["ID"][0], "", "")


def glimmer_predictions(in_handle, ref_recs):
    """Parse Glimmer output, generating SeqRecord and SeqFeatures for predictions
    """
    for rec in GFF.parse(in_handle, target_lines=1000, base_dict=ref_recs):
        yield rec

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit()
    main(*sys.argv[1:])
