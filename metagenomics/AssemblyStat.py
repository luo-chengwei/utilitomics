#! /usr/bin/python

import sys
import os
import re
from Bio import SeqIO

def assemblyStat(infile):
	lengths=[]
	GC_count = 0
	for record in SeqIO.parse(infile, 'fasta'):
		seq=str(record.seq)
		GC_count += seq.upper().count('G')
		GC_count += seq.upper().count('C')
		lengths.append(len(seq))
	sortedLengths = sorted(lengths)
	max_contig = sortedLengths[-1]
	totalLength = sum(lengths)
	n25 = 0
	n50 = 0
	n75 = 0
	acc_length = 0
	for x in sortedLengths:
		acc_length += x	
		if acc_length >= 0.25 * totalLength and n25 == 0:
			n25 = x
		if acc_length >= 0.5 * totalLength and n50 == 0:
			n50 = x
		if acc_length >= 0.75 * totalLength and n75 == 0:
			n75 = x
	GC = 100 * float(GC_count)/totalLength
	return (GC, totalLength, len(lengths), lengths[0], max_contig, n25, n50, n75)


