#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'Chengwei Luo (cluo@broadinstitute.org)'
__version__ = '0.1.0'
__date__ = 'Oct 2013'

"""
getAllReads takes a mapping BAM file and tranform it into *.PE.fa

for usage, type:
python getAllReads.py --help
"""

USAGE = \
"""
getAllReads takes a mapping BAM file and tranform it into *.PE.fa

for usage, type:
python getAllReads.py --help
"""

import sys
import os
import glob
from optparse import OptionParser, OptionGroup


def main(argv = sys.argv[1:]):
	parser = OptionParser(usage = USAGE, version = "Version: " + __version__)
	
	reqOpts = OptionGroup(parser, "Required Options", \
				"These options are required, you may supply them in any order.")
	
	reqOpts.add_option('-i', '--infile', metavar = "FILE",
				help = "The input BAM mapping file [sorted by reads].")
				
	reqOpts.add_option('-o', '--out_prefix', metavar = "STRING",
				help = "The outputfile prefix")
	
	reqOpts.add_option('-f', '--format', metavar = "STRING", default = 'fastq',
				help = "The output file format, has to be \"fastq\", \"fasta\", or \"interleaved_fasta\"")
	
				
	parser.add_option_group(reqOpts)
	
	(options, args) = parser.parse_args(argv)
	
	if options.format not in ['fastq', 'fasta', 'interleaved_fasta']:
		parser.error('[FATAL]: the file format has to be \"fastq\", \"fasta\", or \"interleaved_fasta\"')
		exit(1)

	if options.infile is None:
		parser.error("Infile must be supplied!")
		exit(1)
	
	if options.out_prefix is None:
		parser.error("Outfile prefix must be supplied!")
		exit(1)
		
	bamFile = options.infile
	bamfh = pysam.Samfile(bamFile,'rb')
	for read in readIter:
		if read.is_read1:
			sys.stdout.write('>%s\n%s\n' % (read.qname, read.seq))
		elif read.is_read2:
			sys.stdout.write('>%s\n%s\n' % (read.qname, read.seq))
	ofh1.close()
	ofh2.close()

	bamfh.close()

	
if __name__ == '__main__':
	main()
			