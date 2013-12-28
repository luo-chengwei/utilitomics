#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'Chengwei Luo (cluo@broadinstitute.org)'
__version__ = '0.1.0'
__date__ = 'Dec 2013'


USAGE = \
"""
refMapping2readBins.py takes a mapping against reference genomes BAM file
and split it down to bins of reads belong to certain level of taxa (e.g., all species)

for usage, type:
python refMapping2readBins.py --help
"""

import sys
import os
import re
import glob
from optparse import OptionParser, OptionGroup
import pysam
import taxonomy

parser = OptionParser(usage = USAGE, version = 'Version: ' + __version__)

reqOpts = OptionGroup(parser, "Required Options",
				"These options are required, you may supply them in any order.")

reqOpts.add_option('-b', '--bam', metavar = 'FILE', 
				help = 'input BAM file, sorted and indexed.')
reqOpts.add_option('-o', '--out_dir', metavar = 'DIR',
				help = 'output directory where all the read bins will be created.')
reqOpts.add_option('-t', '--tax_dir', metavar = 'DIR',
				help = 'the directory of the taxonomy dump databases.')
				
parser.add_option_group(reqOpts)

(options, args) = parser.parse_args(sys.argv[1:])

def main():
	if options.bam == None or options.out_dir == None or options.tax_dir == None:
		sys.stderr.write('FATAL: argument input error.\n')
		exit(0)
	
	# initialize taxonomy path
	sys.stdout.write('Initializing taxonomy tree...\n')
	tTree = taxonomy.taxonTree(options.tax_dir)
	sys.stdout.write('Done.\n')
	
	# init output space
	sys.stdout.write('Initializing output dir...\n')
	if not os.path.exists(options.out_dir):
		os.mkdir(options.out_dir)
	else:
		import shutil
		shutil.rmtree(options.out_dir)
		os.mkdir(options.out_dir)
	sys.stdout.write('Done.\n')
	
	
	sys.stdout.write('Initializing taxonomy dict...\n')
	taxonomy_dict = {}
	bfh = pysam.Samfile(options.bam, 'rb')
	for ref in bfh.references:
		seq_id, taxid, type = ref.split('.')
		if type != 'g':
			continue
		if taxid not in taxonomy_dict:
			species_id, species_name = tTree.getRankWithTaxonID(taxid, 'species')
			taxonomy_dict[taxid] = [species_id, species_name]
	bfh.close()
	sys.stdout.write('Done.\n')
		
		
	sys.stdout.write('Categorizing reads...\n')	
	ofhs = {}
	ofhs['plasmid'] = open(options.out_dir + '/plasmids.txt', 'a')
	ofhs['virus'] = open(options.out_dir + '/virus.txt', 'a')
	ofhs['unmapped'] = open(options.out_dir + '/unmapped.txt', 'a')
	ofhs['unknown'] = open(options.out_dir + '/unknown.txt', 'a')
	ofhs['species'] = {}
	
	
	bfh = pysam.Samfile(options.bam, 'rb')
	for read in bfh.fetch(until_eof = True):
		if read.is_unmapped:
			if read.is_read1:
				ofhs['unmapped'].write('%s#1\t%s\t%s\n' % (read.qname, read.seq, read.qual))
			else:
				ofhs['unmapped'].write('%s#2\t%s\t%s\n' % (read.qname, read.seq, read.qual))
			continue
			
		if read.mapq < 30: continue
		
		seq_id, taxid, type = bfh.getrname(read.tid).split('.')
		
		if type == 'p':
			if read.is_read1:
				ofhs['plasmid'].write('%s#1_%s\t%s\t%s\n' % (read.qname, taxid, read.seq, read.qual))
			else:
				ofhs['plasmid'].write('%s#2_%s\t%s\t%s\n' % (read.qname, taxid, read.seq, read.qual))
			continue
			
		if type == 'v':
			if read.is_read1:
				ofhs['virus'].write('%s#1_%s\t%s\t%s\n' % (read.qname, taxid, read.seq, read.qual))
			else:
				ofhs['virus'].write('%s#2_%s\t%s\t%s\n' % (read.qname, taxid, read.seq, read.qual))
			continue
		
		species_id, species_name = taxonomy_dict[taxid]
		if species_id == None:
			if read.is_read1:
				ofhs['unknown'].write('%s#1_%s\t%s\t%s\n' % (read.qname, taxid, read.seq, read.qual))
			else:
				ofhs['unknown'].write('%s#2_%s\t%s\t%s\n' % (read.qname, taxid, read.seq, read.qual))
			continue
			
		outfile = options.out_dir + '/species.' + str(species_id) + '.' + species_name.replace(' ', '_', 100) + '.txt'
		if species_id not in ofhs['species']:
			try:
				ofhs['species'][species_id] = open(outfile, 'a')
			except OSError:
				for sp in ofhs['species']:
					if sp == species_id:
						continue
					if ofhs['species'][sp] != None:
						ofhs['species'][sp].close()
						ofhs['species'][sp] = None			
						ofhs['species'][species_id] = open(outfile, 'a')
						break
		elif ofhs['species'][species_id] == None:
			for sp in ofhs['species']:
				if sp == species_id:
					continue
				if ofhs['species'][sp] != None:
					ofhs['species'][sp].close()
					ofhs['species'][sp] = None			
					ofhs['species'][species_id] = open(outfile, 'a')
					break
		ofh = ofhs['species'][species_id]
		
		if read.is_read1:
			ofh.write('%s#1_%s\t%s\t%s\n' % (read.qname, taxid, read.seq, read.qual))
		else:
			ofh.write('%s#2_%s\t%s\t%s\n' % (read.qname, taxid, read.seq, read.qual))
		
	bfh.close()
	sys.stdout.write('Done.\n')
		
if __name__ == '__main__':
	main()