import sys

DictCod = {}
inF = open('nsl_p1.txt', "r")
for line in inF:
        if line[0] =="#":
                continue

        raw = line.strip().split(":")
        DictCod[raw[0]] = float(raw[1])
inF.close()


a = open(sys.argv[1], "r")
sumN = 0.0
sumA = 0.0
for l in a:
        if l[0] == "#":
                continue
        r = l.split()
        r[3] = r[3].upper()
        if not r[3] in DictCod:
                #print(l) -> very few "Ns" in the dataset
                continue
        sumN += DictCod[r[3]]
        sumA += 3.0
sumS = sumA - sumN
print(str(sumN) + '\t' + str(sumS) + '\t' + str(sumA))
a.close()
