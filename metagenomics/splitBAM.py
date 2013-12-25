#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'Chengwei Luo (cluo@broadinstitute.org)'
__version__ = '0.1.0'
__date__ = 'Oct 2013'

"""
splitBAM takes a mapping BAM file and tranform it into unmapped.PE.[1|2].fa and unmapped.SE.fa

for usage, type:
python getUnmappedReads.py --help
"""

USAGE = \
"""
Usage: %prog [options]

splitBAM takes a mapping BAM file and splits it into designated number of fastq files

for more details, type:
python getUnmappedReads.py --help
"""

import sys
import os
import glob
from optparse import OptionParser, OptionGroup
import pysam


def main(argv = sys.argv[1:]):
	parser = OptionParser(usage = USAGE, version = "Version: " + __version__)
	
	reqOpts = OptionGroup(parser, "Required Options", \
				"These options are required, you may supply them in any order.")
	
	reqOpts.add_option('-i', '--infile', metavar = "FILE",
				help = "The input BAM mapping file [sorted].")
				
	reqOpts.add_option('-o', '--output_prefix', metavar = "STRING",
				help = "The prefix of output files.")
				
	reqOpts.add_option('-n', '--number_chunks', metavar = "INT", type = "int", 
				help = "Number of chunks you want to split the BAM into.")
				
	reqOpts.add_option('-s', '--samtools', metavar = "STRING",
				help = "Location of samtools.")
				
	parser.add_option_group(reqOpts)
	
	(options, args) = parser.parse_args(argv)
	
	# change to interleaved fastq first
	tempfile = options.output_prefix + '.temp'
	os.system('%s view %s | cut -f 1,10,11 > %s' % (options.samtools, options.infile, tempfile))
	
	# output to ofhs
	ofhs = []
	for i in range(options.number_chunks):
		ofhA = open(options.output_prefix + '.' + str(i+1) + '.1.fastq', 'a')
		ofhB = open(options.output_prefix + '.' + str(i+1) + '.2.fastq', 'a')
		ofhs.append((ofhA, ofhB))
		
	tfh = open(tempfile, 'r')
	index = 0
	while 1:
		line = tfh.readline().rstrip('\n')
		if not line:
			break
		tag1, seq1, qual1 = line.split('\t')
		line = tfh.readline().rstrip('\n')
		tag2, seq2, qual2 = line.split('\t')
		i = index % options.number_chunks
		index += 1
		ofhs[i][0].write('@%s\n%s\n+\n%s\n' % (tag1, seq1, qual1))
		ofhs[i][1].write('@%s\n%s\n+\n%s\n' % (tag2, seq2, qual2))
	tfh.close()
	
	os.remove(tempfile)
	
	for ofhA, ofhB in ofhs:
		ofhA.close()
		ofhB.close()
	
if __name__ == '__main__':
	main()