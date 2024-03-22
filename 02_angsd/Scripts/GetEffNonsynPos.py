import sys
import itertools
from Bio import SeqIO, Seq

def getCodonNsynPos():
        DictTr = {}
        Nucl = ["A", "C", "T", "G"]
        Cods = itertools.product("ACTG", repeat=3)
        for Cod in Cods:
                C = Cod[0] + Cod[1] + Cod[2]
                fCodSeq = Seq.Seq(C)
                P = fCodSeq.translate()
                DictTr[C] = P


        DictCod = {}
        Cods = itertools.product("ACTG", repeat=3)
        for Cod in Cods:
                C = Cod[0] + Cod[1] + Cod[2]
                P = Seq.Seq(C).translate()

                n1 = 0
                n2 = 0
                n3 = 0
                for n in Nucl:
                        if not n == Cod[0]:
                                C1 = n +  Cod[1] + Cod[2]
                                P1 = Seq.Seq(C1).translate()
                                if not DictTr[C] == DictTr[C1]:
                                        n1 += 1
                        if not n == Cod[1]:
                                C2 = Cod[0] + n + Cod[2]
                                P2 = Seq.Seq(C2).translate()
                                if not DictTr[C] == DictTr[C2]:
                                        n2 += 1

                        if not n == Cod[2]:
                                C3 = Cod[0] + Cod[1] + n
                                P3 = Seq.Seq(C3).translate()
                                if not DictTr[C] == DictTr[C3]:
                                        n3 += 1
                DictCod["".join(Cod)] = (n1, n2, n3)

        return DictCod




##print(getCodonNsynPos())
DictCod = getCodonNsynPos()
a = open(sys.argv[1], "r")

for line in a:
        if line[0] == "#":
                continue

        raw = line.strip().split()
        C = raw[3].upper()
        if not C in DictCod:
                continue
        Sites = DictCod[C]
        for i in [0,1,2]:
                print(raw[i] + "\t" + str( round(float(Sites[i])/3,2) ))


a.close()
