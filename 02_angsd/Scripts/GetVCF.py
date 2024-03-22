##The script assumes the same scaffold for both files
import sys

Cods = open(sys.argv[1], "r")
vcf = open(sys.argv[2], "r")

DictP = {}
for line in Cods:
        if line[0] == "#":
                continue

        raw = line.strip().split()
        #print(raw)
        if len(raw) < 3:
                continue
        DictP[str(int(raw[0])+1)] = ""
        DictP[str(int(raw[1])+1)] = ""
        DictP[str(int(raw[2])+1)] = ""

Cods.close()

for l in vcf:
        if l[0] == "#":
                print(l.strip())
                continue
        raw = l.split()
        #print(raw)
        if raw[1] in DictP:
                print(l.strip())
vcf.close()
