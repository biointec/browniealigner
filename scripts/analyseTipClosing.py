import sys
import re
import hashlib
from itertools import izip
from Bio import SeqIO
maxReadLength=300
class SAMEntry:
    def __init__(self, props):
	self.qname = props[0]
	self.flag = props[1]
	self.rname = props[2]
	self.pos = int(props[3]) - 1 #change to 0-based
	self.mapq = props[4]
	self.cigar = re.findall(r'(\d+)([A-Z]{1})', props[5])
	self.cigar = [(int(c), s) for c, s in self.cigar]
	self.rnext = props[6]
	self.pnext = props[7]
	self.tlen = props[8]
	self.seq = props[9]
	self.qual = props[10]
	self.optional = {}
	for i in range(11, len(props)):
	    field = props[i].split(':')
	    if field[1] == 'i':
		field[2] = int(field[2])
	    elif field[1] == 'f':
		field[2] = float(field[2])
	    self.optional[field[0]] = [field[1], field[2]]
	self.len = len(self.seq)
	if self.cigar[0][1] == 'S':
	    self.len -= self.cigar[0][0]
	if self.cigar[-1][1] == 'S':
	    self.len -= self.cigar[-1][0]
def sam_parse(iFName):
    infile = open(iFName)
    sam = []
    for line in infile:
	if line.startswith('@'): #header
	    continue
	else: #entry
	    props = line.strip().split('\t')
	    if props[5] != '*':
		sam += [SAMEntry(props)]
    infile.close()
    return sam
def defineStatistics(samRecords):
    i = 0;
    TP =0
    FP =0 
    TN=0
    FN =0
    while (i < len(samRecords)-1):
	while (not samRecords[i].qname.endswith("1")):
	    i = i +1    
	first = samRecords[i]
	while (not samRecords[i].qname.endswith("2")):
	    i = i +1    
	second = samRecords[i]
	
	if (abs (second.pos - first.pos) >500) :
	    if (first.qname.split("_")[1]=='0'):
		TN = TN +1
		#print ("************** TN **************")
                #print (abs (second.pos - first.pos))
                #print (">"+first.qname)
                #print (first.seq)
                #print (">"+second.qname)
                #print (second.seq)

	    else:
		FP = FP +1
		#print ("************** FP **************")
		#print (abs (second.pos - first.pos))
		#print (">"+first.qname)
		#print (first.seq)
		#print (">"+second.qname)
		#print (second.seq)
		
	else:
	    if (samRecords[i].qname.split("_")[1]=='0'):
		#print ("************** FN **************")
		#print (abs (second.pos - first.pos))
		#print (">"+second.qname)
		#print (first.seq)
		#print (">"+second.qname)
		#print (second.seq)
		FN = FN +1
	    else:
		TP = TP +1

	i = i +2
    print ("Number of True Positive :" +str (TP))
    print ("Number of False Positive :" +str (FP))
    print ("Number of True Negative :" +str (TN))
    print ("Number of False Negative :" +str (FN))
    print ("Gain :" + str( (TP-FP)*100/(TP+FN)) +"%")
    print ("Sensitivity :" + str( (TP)*100/(TP+FN)) +"%")
    print ("FP Rate :" + str( (FP)*100/(TN+FP)) +"%")
def main(argv=None):
    if argv == None:
	argv = sys.argv
    if len(argv) < 1:
	print('Usage: analyseTipClosing.py samFile.sam')
        exit()       
    samFile = argv[1]
    records = sam_parse(samFile)
    defineStatistics(records)
    
main()
