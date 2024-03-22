import sys

a = open(sys.argv[1], "r")
for line in a:
        if line[0] == "#":
                print(line.strip())
                continue

        raw = line.split("\t")
        n = raw[0].split(":")[1].split("|")[0]
        raw[0] = n
        print("\t".join(raw))
a.close()
