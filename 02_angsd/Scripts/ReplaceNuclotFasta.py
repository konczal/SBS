from Bio import SeqIO
from Bio.Seq import Seq

MajorF = open("SBS_freqs.mafs", "r")

DictMaj = {}
for l in MajorF:
        raw = l.strip().split()
        if raw[0] == "chromo":
                continue
        raw = l.strip().split()
        if not raw[2] == raw[4]:
                DictMaj[(raw[0], raw[1])] = l

MajorF.close()

DictFasta = {}
DictList = {}
data = "/media/raid/home/mkonczal/Projects/SBS_2021/00_data/genome/SBS_final.scaffolds.fasta"
seqOrder = []
for record in SeqIO.parse(data, "fasta"):
        seqOrder.append(record.id)
        DictFasta[record.id]  = record
        DictList[record.id] = list(str(record.seq))
#counter = 1
for k in DictMaj:
        pos = int(k[1]) - 1
        Maj = DictMaj[k].split()[2]
        if not DictList[k[0]][pos] == Maj:
                DictList[k[0]][pos] = Maj
                #print(str(counter))
                #counter += 1

for s in DictFasta:
        DictFasta[s].seq = Seq("".join(DictList[s]))

with open("SBS_final.scaffolds_SBSmajor.fasta", "w") as output_handle:
        for key in seqOrder:
                SeqIO.write(DictFasta[key], output_handle, "fasta")
