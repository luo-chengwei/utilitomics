import sys
import os
import re

infile1 = sys.argv[1]
infile2 = sys.argv[2]
readList = sys.argv[3]
outfile1 = sys.argv[4]
outfile2 = sys.argv[5]

blackList = {}
for read in open(readList,'r'):
	blackList[read.rstrip('\n')] = 1

ifh1 = open(infile1, 'r')
ifh2 = open(infile2, 'r')
ofh1 = open(outfile1, 'w')
ofh2 = open(outfile2, 'w')

while 1:
	tag1 = ifh1.readline().rstrip('\n')
	if not tag1:
		break
	seq1 = ifh1.readline().rstrip('\n')
	dis1 = ifh1.readline().rstrip('\n')
	qual1 = ifh1.readline().rstrip('\n')

	tag2 = ifh2.readline().rstrip('\n')
	seq2 = ifh2.readline().rstrip('\n')
	dis2 = ifh2.readline().rstrip('\n')
	qual2 = ifh2.readline().rstrip('\n')

	if tag1[1:] in blackList:
		continue

	ofh1.write('%s\n%s\n%s\n%s\n' % (tag1, seq1, dis1, qual1))
	ofh2.write('%s\n%s\n%s\n%s\n' % (tag2, seq2, dis2, qual2))

ifh1.close()
ifh2.close()
ofh1.close()
ofh2.close()

