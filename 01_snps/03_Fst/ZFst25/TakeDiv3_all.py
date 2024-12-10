'''
reated 27.01.2017

@author: mkonczal
'''

import sys
import fileinput

HOM = ["0/0","1/1","2/2","3/3", "4/4"]

SBS = [9,10,11,12,13,14,15]
RNS = [16,17,18,19,20,21,22,23,24]
lSBS_h = 2*(len(SBS))
lRNS_h = 2*(len(RNS))

CMIN = 25
CSUB = 26
CCAN = 27

#[polymorphism, substiutions] for each class of SNPs
EFFECT_ORDERED = ["stop_gained", "stop_lost", "start_lost", "missense_variant", "synonymous_variant", "intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "intragenic_variant", "downstream_gene_variant", "upstream_gene_variant", "intergenic_region"]
EFFECTS_SBS = {"stop_lost":[0,0], "start_lost":[0,0], "intragenic_variant":[0,0], "stop_gained":[0,0], "3_prime_UTR_variant":[0,0], "5_prime_UTR_variant":[0,0], "downstream_gene_variant":[0,0], "intergenic_region":[0,0], "intron_variant":[0,0],"missense_variant":[0,0],"splice_region_variant":[0,0],"synonymous_variant":[0,0], "upstream_gene_variant":[0,0]}
EFFECTS_RNS = {"stop_lost":[0,0], "start_lost":[0,0], "intragenic_variant":[0,0], "stop_gained":[0,0], "3_prime_UTR_variant":[0,0], "5_prime_UTR_variant":[0,0], "downstream_gene_variant":[0,0], "intergenic_region":[0,0], "intron_variant":[0,0],"missense_variant":[0,0],"splice_region_variant":[0,0],"synonymous_variant":[0,0], "upstream_gene_variant":[0,0]}

OTHER_SBS = [0,0]
OTHER_RNS = [0,0]

POLI_SBSRNS = 0
######
comSNP_SBS=0   ##int(sys.argv[1])   2 #<-----Should be changed....
comSNP_RNS=0   ##int(sys.argv[2])   2
######


nSNPs = 0
nSNPsdiv = 0
nSNPspoly = 0
SBSsub = 0
RNSsub = 0

for line in fileinput.input():
	if line[0] == "#":
		continue
	row = line.split("\t")


####FILTER INDELS#########
#	if not len(row[3]) == 1:
#		continue
# 
#	if not len(row[4].split(",")[0]) == 1:
#		continue
#######


##!!27.01.2017 <- take only one class, the highest in the "rank"
####DEFINE SNP CLASS###########
	SNP_CLASS = []
	for keyEFF in EFFECT_ORDERED:
		if  keyEFF in row[7]:
			SNP_CLASS.append(keyEFF)
			break
##!!
########TAKE ANCESTRAL##########
	if len(SNP_CLASS) > 1:
		continue

	out = row[CCAN].split(":")[0]
	if not out.split("/")[0] == out.split("/")[1]:
		continue
	ANC = out.split("/")[0]

#########TAKE FREQS FOR POPS####
	ANCcount_SBS = 0
	ANCcount_RNS = 0

	for iSBS in SBS:
		ANCcount_SBS += row[iSBS].split(":")[0].count(ANC)

	for iRNS in RNS:
		ANCcount_RNS += row[iRNS].split(":")[0].count(ANC)


	#print str(ANCcount_RNS)


##!!27.01
#########EXCLUDE POLIMORPHIC IN TWO SPECIES#########
	if min(ANCcount_SBS, lSBS_h - ANCcount_SBS) > 0:
		if  min(ANCcount_RNS, lRNS_h - ANCcount_RNS) > 0: #<----polimorphic in both species
			POLI_SBSRNS += 1
			continue
	

#########EXCLUDE IF SITES DIFFER BETWEEN NOT POLIMORPHIC SPECIES#############
	if min(ANCcount_SBS, lSBS_h - ANCcount_SBS) > 0:
		if not ANCcount_RNS == lRNS_h:
			continue

	if min(ANCcount_RNS, lRNS_h - ANCcount_RNS) > 0:
		if not ANCcount_SBS == lSBS_h:
			continue		



#######DEFINE EFFECT###########

	POL = [0,0]
	if (comSNP_SBS < lSBS_h - ANCcount_SBS and min(ANCcount_SBS, lSBS_h - ANCcount_SBS) > 0):
		POL[0] += 1

	if (comSNP_RNS < lRNS_h - ANCcount_RNS and min(ANCcount_RNS, lRNS_h - ANCcount_RNS) > 0):
		POL[1] += 1
		
	if SNP_CLASS == []:
		OTHER_SBS[0] += POL[0]
		OTHER_RNS[0] += POL[1]
			
	if not SNP_CLASS == []:
		for keyEFF_B in SNP_CLASS:
			EFFECTS_SBS[keyEFF_B][0] += POL[0]
			EFFECTS_RNS[keyEFF_B][0] += POL[1]

	nSNPs += 1
	
#####IF POLYMORPHIC IN AT LEAST ONE SPECIES: CONTINUE
	if not min(ANCcount_SBS, lSBS_h-ANCcount_SBS) == 0:
		nSNPspoly += 1
		continue
	if not min(ANCcount_RNS, lRNS_h-ANCcount_RNS) == 0:
		nSNPspoly += 1
		continue


###DEFINE DIVERGENCE
	##print str(ANCcount_SBS) + "\t" + str(ANCcount_RNS)
	if (ANCcount_SBS == 0 and ANCcount_RNS == lRNS_h):  ###<---SBS substitution 
		if SNP_CLASS == []:
			OTHER_SBS[1] += 1
		else:
			for keyEFF_D in SNP_CLASS:
				EFFECTS_SBS[keyEFF_D][1] += 1
			
		SBSsub += 1
		nSNPsdiv += 1
		continue
		
	if (ANCcount_SBS == lSBS_h and ANCcount_RNS == 0): ###<---RNS substitution
		if SNP_CLASS == []:
			OTHER_RNS[1] += 1

		else:
			for keyEFF_F in SNP_CLASS:
				EFFECTS_RNS[keyEFF_F][1] += 1
				 
		nSNPsdiv += 1
		RNSsub += 1

#print "Effect" + "\tPoly_SBS\tSubs_SBS\tPoly_RNS\tSubs_RNS"
#for keyEFF_C in EFFECTS_SBS:
#	print keyEFF_C + "\t" + "\t".join(map(str, EFFECTS_SBS[keyEFF_C])) + "\t" + "\t".join(map(str, EFFECTS_RNS[keyEFF_C])) 
#print "other_effects" + "\t" + "\t".join(map(str, OTHER_SBS)) + "\t" + "\t".join(map(str, OTHER_RNS))
#print "\n#####OTHER COUNTS####\npolimorphic_both\t" + str(POLI_SBSRNS)
#print "Number of all SNPs\t" + str(nSNPs)
#print "Number of SNPs for divergence\t" + str(nSNPsdiv)
#print "Number of SNPs polimorphic in one species\t" + str(nSNPspoly)
print str(SBSsub) + "\t" + str(RNSsub)
