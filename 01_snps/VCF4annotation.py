'''
Created on 18 de may. de 2016

@author: mkonczal
'''
import sys
import fileinput

def FilterVCF():
    #a = open(infile, "r")
    for l in fileinput.input():
        line = l.strip()
        if line[0] == "#":
            if line[:9] == "##contig=":
                print line.split("|")[0] + "," + ">"
                continue
            print line

        else:
            row = line.split("\t")
            sc = row[0].split("|")[0]
            print sc + "\t" + "\t".join(row[1:])

if __name__ == '__main__':
    #FilterVCF(sys.argv[1])
    FilterVCF()

