import sys

cod = open(sys.argv[1], "r")
vcf = open(sys.argv[2], "r")

DictPos = {}
for line in vcf:
        if line[0] == "#":
                continue

        raw = line.split()
        DictPos[int(raw[1])] = ""
vcf.close()

for l in cod:
        if l[0] == "#":
                continue
        r =l.split()
        if (int(r[0])+1 in DictPos) and (int(r[1])+1 in DictPos) and (int(r[2])+1 in DictPos):
                print(l.strip())
cod.close()
