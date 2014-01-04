#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'Chengwei Luo (cluo@broadinstitute.org)'
__version__ = '0.1.0'
__date__ = 'Dec 2013'


USAGE = \
"""
dynamic_assembler.py takes the directory from refMapping2readBins.py, the ref. database fasta
file, and the taxonomy name dmp file (NCBI ftp://ftp.ncbi.nlm.nih.gov:/pub/taxonomy/taxonomy.tar.gz)

It requires velvet and spades assemblers being in $ENV

The output of this script would be assembled contigs.

for usage, type:
python dynamic_assembler.py --help
"""

import sys
import os
import re
import glob
from optparse import OptionParser, OptionGroup
import pysam
from Bio import SeqIO
import taxonomy

parser = OptionParser(usage = USAGE, version = 'Version: ' + __version__)

reqOpts = OptionGroup(parser, "Required Options",
				"These options are required, you may supply them in any order.")

reqOpts.add_option('-i', '--in_dir', metavar = 'DIR', 
				help = 'The directory output by refMapping2readBins.py.')
reqOpts.add_option('-r', '--ref', metavar = 'FILE',
				help = 'The reference database used by ref. mapping')
reqOpts.add_option('-t', '--tax_dir', metavar = 'DIR',
				help = 'the directory of the taxonomy dump databases.')
reqOpts.add_option('-o', '--out_dir', metavar = 'DIR',
				help = 'the output directory where all the output files could be found.')		
	
parser.add_option_group(reqOpts)

(options, args) = parser.parse_args(sys.argv[1:])

def main():
	if options.in_dir == None or options.out_dir == None \
			or options.tax_dir == None or options.ref == None:
		sys.stderr.write('FATAL: argument input error.\n')
		exit(0)
		
	
		
if __name__ == '__main__':
	main()
		
	