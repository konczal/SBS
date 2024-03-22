import sys

sfs = open(sys.argv[1], "r")
ID = sys.argv[2]

line =sfs.readline()
raw = line.strip().split()
if not len(raw) == 3:
        print("ERROR " + ID + "\t" + line)

HET = float(raw[1])/(float(raw[0])  + float(raw[1]) + float(raw[2]))
print(ID + "\t" + str(round(HET,6)) + "\t" + str(round(float(raw[1]),2)))
sfs.close()
