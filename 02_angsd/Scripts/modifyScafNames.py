import sys

DictSc = {}

a = open("Scaffold_names_dict.txt", "r")
for line in a:
        raw = line.strip().split()
        DictSc[raw[1]] = raw[0]
a.close()

b = open(sys.argv[1], "r")
for l in b:
        r = l.strip().split()
        nr = DictSc[r[0]]
        r[0] = nr
        nnr = str(int(r[1]) + 1)
        r[1] = nnr
        print("\t".join(r))
b.close()
