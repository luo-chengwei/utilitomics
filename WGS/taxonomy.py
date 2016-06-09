#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This library works with the taxonomy db download from NCBI FTP:

ftp://ftp.ncbi.nih.gov:/pub/taxonomy/taxdump.tar.gz

basically you can init an taxonomy tree obj by:

tTree = taxonomy.taxonTree(dir)   # tree being the tax dump db unzip dir
then you can basically get the taxonomy ranks using:

path = tTree.getNamePathWithTaxonID(taxonID)

"""

import sys
import glob
import networkx as nx

class taxonTree:
	def __init__(self):
		self.tree = nx.DiGraph()
		self.name2taxonID = {}
		self.taxonID2name = {}
	
	def __init__(self, dir):
		self.tree = nx.DiGraph()
		try:
			nameslib = glob.glob(dir + '/names.dmp')[0]
		except:
			sys.stderr.write('FATAL: cannot find the names.dmp file, initilization aborted.\n')
			exit(1)
		
		try:
			nodeslib = glob.glob(dir + '/nodes.dmp')[0]
		except:
			sys.stderr.write('FATAL: cannot find the nodes.dmp file, initilization aborted.\n')
			exit(1)
		
		# construct name -> taxonID mapping	
		self.name2taxonID = {}
		self.taxonID2name = {}
		for line in open(nameslib, 'r'):
			cols = line.split('\t')
			taxonID = cols[0]
			name = cols[2]
			self.name2taxonID[name] = taxonID
			if cols[-2] == 'scientific name':
				self.taxonID2name[taxonID] = name
			
		# construct node tree
		edges = []
		nodes = {}
		for line in open(nodeslib, 'r'):
			cols = line.split('\t')
			parentNode = cols[2]
			childrenNode = cols[0]
			rank = cols[4]
			nodes[childrenNode] =rank
			if childrenNode != parentNode:
				edges.append((childrenNode, parentNode))
		
		self.tree.add_edges_from(edges)
		nx.set_node_attributes(self.tree, 'rank', nodes)
		
	def getTaxonIDPathWithTaxonID(self, taxonID):
		path = [taxonID]
		currentNode = taxonID
		if currentNode not in self.tree.nodes():
			return []
		while len(self.tree.successors(currentNode)) > 0:
			path.append(self.tree.successors(currentNode)[0])
			currentNode = self.tree.successors(currentNode)[0]
		return path
		
	def getTaxonIDPathWithName(self, name):
		if name not in self.name2taxonID:
			return []
		path = self.getTaxonIDPathWithTaxonID(self.name2taxonID[name])
		return path
		
	def getNamePathWithTaxonID(self, taxonID):
		taxonIDPath = self.getTaxonIDPathWithTaxonID(taxonID)
		namePath = []
		for x in taxonIDPath:
			rank = self.tree.node[x]['rank']
			try:
				name = self.taxonID2name[x]
			except KeyError:
				name = ''
			namePath.append((name, rank))
		return namePath
		
	def getNamePathWithName(self, name):
		if name not in self.name2taxonID:
			return []
		path = self.getNamePathWithTaxonID(self.name2taxonID[name])
		return path
	
	def getRankWithTaxonID(self, taxonID, rank):
		taxonIDPath = self.getTaxonIDPathWithTaxonID(taxonID)
		namePath = self.getNamePathWithTaxonID(taxonID)
		for tid, (name, x) in zip(taxonIDPath, namePath):
			if x == rank:
				return tid, name
		return None, None
				
	def getRankWithName(self, name, rank):
		if name not in self.name2taxonID:
			return None, None
		return self.getRankWithTaxonID(self.name2taxonID[name], rank)
		
# End of class