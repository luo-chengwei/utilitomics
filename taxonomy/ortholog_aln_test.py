#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'Chengwei Luo (cluo@broadinstitute.org)'
__version__ = '0.1.0'
__date__ = 'December 2013'


USAGE = \
"""
This script calculates the significance of orthologs that differentiate categories.

This script takes a multi-newick tree file, in which each leaf is a gene sequence categorized by
the sample, the category, and the ortholog group it belongs to; ecach entry should
looks like this:

gene_id|sample_id|category_id|ortholog_id

Generally, you use it like this:

$python ortholog_aln_test -i [input_aln_file] -o [output_file] [options]

Dependencies:

- Biopython

for detailed usage, type:
python ortholog_aln_test.py --help
"""

import sys
import os
import re
import glob
from optparse import OptionParser, OptionGroup
import multiprocessing as mp
from Bio import Phylo, SeqIO
from itertools import combinations, permutations
from subprocess import PIPE, Popen, call
from time import ctime, time
from datetime import timedelta
import shutil
from cStringIO import StringIO
from scipy.stats import binom

def quartet_analysis(tree, quartets, perms):
	# generate p-value distribution using Bernoulli model
	numTopo1 = 0
	numTopo2 = 0
	numTopo3 = 0
	for perm in perms:
		nodeA = perm[0][0]
		nodeB = perm[0][1]
		nodeC = perm[1][0]
		nodeD = perm[1][1]
		distAB=tree.distance(nodeA, nodeB)
		distAC=tree.distance(nodeA, nodeC)
		distAD=tree.distance(nodeA, nodeD)
		if distAB == min(distAB, distAC, distAD):
			numTopo1+=1
		elif distAC == min(distAB, distAC, distAD):
			numTopo2+=1
		else:
			numTopo3+=1
	
	P = float(numTopo1) / len(perms)  # P is the binomial dist's p, estimated by permutations
	
	# count diff types of quartets
	if len(quartets) > 1000:
		selQ = random.sample(quartets, 1000)
	else:
		selQ = quartets
		
	numTopo1 = 0
	numTopo2 = 0
	numTopo3 = 0
	for quartet in selQ:
		nodeA = quartet[0][0]
		nodeB = quartet[0][1]
		nodeC = quartet[1][0]
		nodeD = quartet[1][1]
		distAB=tree.distance(nodeA, nodeB)
		distAC=tree.distance(nodeA, nodeC)
		distAD=tree.distance(nodeA, nodeD)
		if distAB == min(distAB, distAC, distAD):
			numTopo1+=1
		elif distAC == min(distAB, distAC, distAD):
			numTopo2+=1
		else:
			numTopo3+=1
	
	# calculate the p using binomial dist.
	p = binom.sf(numTopo1, len(selQ), P)  # one tailed p-val, binomial dist. survival function
	return p

def quartet_test(tree_string):
	tree = Phylo.read(StringIO(tree_string), "newick")
	nodes = {}
	for leaf in tree.get_terminals():
		if leaf.name == None:
				continue
		category = leaf.name.split('|')[-2]
		if category not in nodes:
			nodes[category] = []
		nodes[category].append(leaf)
	
	results = []
	typeCombos = list(combinations(nodes.keys(), 2))
	for category_1, category_2 in typeCombos:
		if len(nodes[category_1]) < 2 or len(nodes[category_2]) < 2:
			continue
		
		# generate combinations for quartet_test
		combo1 = list(combinations(nodes[category_1], 2))
		combo2 = list(combinations(nodes[category_2], 2))
		quartets = []
		for c1 in combo1:
			for c2 in combo2:
				quartets.append((c1, c2))
		
		# generate combinations for permutation test
		perms = []
		combos = list(permutations(nodes[category_1] + nodes[category_2], 4))[:1000]
		for c in combos:
			perms.append(((c[0], c[1]), (c[2], c[3])))
		
		# calculate the p-value
		p = quartet_analysis(tree, quartets, perms)
		
		# append to results
		count_1 = len(nodes[category_1])
		count_2 = len(nodes[category_2])
		results.append((category_1, count_1, category_2, count_2, p))
	
	return results

def ind_test(x):
	ortholog_id, tree = x
	
	# calculate p-value
	xp = quartet_test(tree)
	results = []
	for category_1, count_1, category_2, count_2, p_val in xp:
		results.append((ortholog_id, category_1, count_1, category_2, count_2, p_val))
	return results
	
def BH_correction(p_vals):
	m = len(p_vals)
	adj_pvals = []
	for k, x in enumerate(sorted(p_vals, key = lambda x: x[-1])):
		px = min(1, x[-1]*(float(m)/(k+1)))
		newx=[a for a in x]
		newx.append(px)
		adj_pvals.append(newx)
	return adj_pvals

	
def partition_list(x, l):
	for i in range(0, len(x), l):
		yield x[i:i+l]

def main(argv = sys.argv[1:]):
	parser = OptionParser(usage = USAGE, version = "Version: " + __version__)
	
	reqOpts = OptionGroup(parser, "Required Options",
				"These options are required, you may supply them in any order.")
	
	reqOpts.add_option('-i', '--infile', type = "string", metavar = "FILE",
				help = "The input aln file.")
				
	reqOpts.add_option('-o', '--outfile', type = "string", metavar = "FILE",
				help = "The outpu tfile")
				
	parser.add_option_group(reqOpts)
	
	optOpts = OptionGroup(parser, "Optional Options", 
				"These options are optional, you may supply them in any order.")
	
	optOpts.add_option('-t', '--num_proc', type = "int", metavar = "INT",
				help = "Number of processors available [default: 1]")
				
	parser.add_option_group(optOpts)
	
	(options, args) = parser.parse_args(argv)
	
	
	if options.infile is None:
		parser.error("Input aln file must be supplied!")
		exit(1)
	
	if options.outfile is None:
		parser.error("Outfile must be supplied!")
		exit(1)
	
	num_proc = options.num_proc
	infile = options.infile
	outfile = options.outfile
	
	## start analysis
	# kickstart
	sys.stdout.write("ortholog_aln_test started at %s\n"%(ctime()))
	sys.stdout.flush()
	
	# partition them into 100 temp parts
	trees = []
	ifh = open(infile, 'r')
	while 1:
		gid = ifh.readline().rstrip('\n')
		if not gid:
			break
		tree = ifh.readline().rstrip('\n')
		trees.append((gid, tree))
	ifh.close()
	
	p_vals = []
	tree_pars = partition_list(trees, 100)
	
	for tree_par in tree_pars:
		cmds = [[gid, tree] for gid, tree in tree_par]
		pool = mp.Pool(num_proc)
		x = pool.map_async(ind_test, cmds)
		pool.close()
		pool.join()
		for r in x.get():
			p_vals += r
	# FDR Correction using Benjamini-Hochberg approach
	q_vals = BH_correction(p_vals)
	
	# output results
	ofh = open(outfile, 'w')
	ofh.write('#ortholog_id\tcategory_1\tcount_1\tcategory_2\tcount_2\tpval\tadj_pval(FDR_BH)\n')
	for q_val in q_vals:
		ortholog_id, category_1, count_1, category_2, count_2, p_val, q_val = q_val
		ofh.write('%s\t%s\t%i\t%s\t%i\t%.4f\t%.4f\n' % \
				(ortholog_id, category_1, count_1, category_2, count_2, p_val, q_val))
	ofh.close()
	
	# end
	sys.stdout.write("ortholog_aln_test finished at %s\n"%(ctime()))
	sys.stdout.flush()
	
if __name__ == '__main__':
	main()
			