#!/usr/bin/env python
# -*- coding: utf-8 -*- 

from __future__ import with_statement

__author__ = 'Chengwei Luo (luo.chengwei@gatech.edu)'
__version__ = '0.1.0'
__date__ = '09-03-2013'

import numpy as np
import sys
import os
import re
import argparse as ap
import subprocess as subp
import tempfile as tf
import matplotlib.pyplot as plt

tax_units = "kpcofgs"
tax_lookup = {'k':'Kingdom', 'p':'Phylum', 'c':'Class', 'o':'Order',
				'f':'Family', 'g':'Genus', 's':'Species'}

def read_params(args):
	p = ap.ArgumentParser( description= 
            "DESCRIPTION\n"
            "MG16S version "+__version__+" ("+__date__+"): Metagenomic paired-end reads 16S portion extractor\n"
            "AUTHORS: "+__author__+"\n\n"
            "COMMON COMMANDS\n\n"
            ".......(under construction)\n\n",
            formatter_class=ap.RawTextHelpFormatter )
    
	arg = p.add_argument
	
	arg( 'in', metavar='INPUT_FILE', type=str, nargs='?', default=None, help= 
         "the input file should be:\n"
         "* a multi-fasta file containing metagenomic reads in this format:\n"
         ">tag1A\nseq1A\n>tag1B\nseq1B\n>tag2A\nseq2A\n>tag2B\n...\n\n"
         "where tag1A and tag1B are the tags of the paired-end reads\n\n" )   
    
	arg( 'outprefix', metavar='OUTPUT_PREFIX', type=str, nargs='?', default=None, help= 
    	"the prefix for all the output files \n")

	arg( '-v','--version', action='store_true', help="Prints the current MG16S version and exit\n" )
	
	arg( '--probe', metavar="", type=str, default = None,
		help = "The BLASTN reference of 16S DNA sequence " )
	
	arg( '--GreenGenesDB', metavar="", type=str, default = None,
		help = "The GreenGene reference database " )
	
	arg( '--GreenGenesTax', metavar="", type=str, default = None,
		help = "The GreenGenes database linking OTU# to taxonomy annotation" )
    
	arg( '--blast', metavar="", type=str, default = None, help =
		'Full path and name of the blastall executable. This option allows \n'
		'MG16S to reach the executable even when it is not in the system \n'
		'PATH or the system PATH is unreachable\n' )
		
	arg( '--formatdb', metavar="", type=str, default = None, help =
		'Full path and name of the formatdb executable. This option allows \n'
		'MG16S to reach the executable even when it is not in the system \n'
		'PATH or the system PATH is unreachable\n' )
		
	arg( '--evalue', metavar="", default="1e-6", type=str, help = 
		"evalue threshold for the blasting\n"
		"[default 1e-6]"   )
    
	arg( '--alignLength', metavar="", default="50", type=str, help =
		"alignment length threshold for the blastn search against GreenGenes\n"
		"[default: 50bp]" )
    	
	arg( '--alignIdentity', metavar="", default="97", type=str, help =
		"alignment nt identity threshold for the blastn search against GreenGenes\n"
		"[default: 97]" )
    	
	arg( '--taxLevel', metavar='TAXONOMIC_LEVEL', type=str, 
		choices=tax_units, default='g', help = 
		"The taxonomic level for the relative abundance output:\n"
		"'k' : kingdoms (Bacteria and Archaea)\n"
		"'p' : phyla\n"
		"'c' : classes\n"
		"'o' : orders\n"
		"'f' : families\n"
		"'g' : genera\n"
		"'s' : species\n"
		"[default 'g']" )
		
	arg( '--nproc', metavar="N", type=int, default=1, help = 
		"The number of CPUs to use for parallelizing the blasting\n"
		"[default 1, i.e. no parallelism]\n" )
    
	return vars(p.parse_args())
## End of read_params

def printOptions(pars):
	sys.stdout.write("\nMG16S version "+__version__+"\t("+__date__+")"+"\n")
	sys.stdout.write("Starting protocol, the options you selected are:\n"
					"=====================================\n")
	sys.stdout.write("--Probe 16S sequence file: "+pars['probe']+"\n"
					"--GreenGenes database file: "+pars['GreenGenesDB']+"\n"
					"--GreenGenes taxonomy file: "+pars['GreenGenesTax']+"\n"
					"--Blast binary: "+(pars['blast'] if pars['blast'] else "default PATH")+"\n"
					"--formatdb binary: "+(pars['formatdb'] if pars['formatdb'] else "default PATH")+"\n")
	sys.stdout.write("--Blast search E-value: "+pars['evalue']+"\n"
					"--Blast alignment length: "+pars['alignLength']+"bp\n"
					"--Blast alignment nt identity: "+pars['alignIdentity']+"%\n")
	sys.stdout.write("--Taxonomy stat level: "+tax_lookup[pars['taxLevel']]+"\n"
					"--Number of processors in use: "+str(pars['nproc'])+"\n"
					"--input file: "+pars['in']+"\n"
					"--output files prefix: "+pars['outprefix']+"\n"
					"=====================================\n")
					
## End of printOptions					

def run_blastn(formatdb_exe, blast_exe, infile, 
				outfile, dbfile, evalue, nproc):
	## formatdb command
	try:
		formatdb_cmd = [ formatdb_exe if formatdb_exe else 'formatdb', 
                        "-pF", "-n"+dbfile, "-i"+dbfile ] 
                        
		retcode = subp.call( formatdb_cmd )
	except OSError:
		sys.stderr.write( "OSError: fatal error running formatdb. Is formatdb in the system path?\n" )
		sys.exit(1)
	except ValueError:
		sys.stderr.write( "ValueError: fatal error running formatdb.\n" )
		sys.exit(1)
	except IOError:
		sys.stderr.write( "IOError: fatal error running formatdb.\n" )
		sys.exit(1)
	if retcode == 13:
		sys.stderr.write( "Permission Denied Error: fatal error running formatdb." 
          "Is the formatdb binary in the path with execution and read permissions?\n" )
		sys.exit(1)
        
    ## blastn command
	try:
#		logfh=open(logfile,'w')
		blast_cmd = [ blast_exe if blast_exe else 'blastall', "-pblastn",
   					"-X150","-q-1","-F F", "-d"+dbfile,"-b1","-v1","-m8","-e"+evalue,
   					"-a"+str(nproc), "-i"+infile, "-o"+outfile]
		retcode = subp.call( blast_cmd, stdout=subp.PIPE)
	except OSError:
		sys.stderr.write( "OSError: fatal error running blastn. Is blastall in the system path?\n" )
		sys.exit(1)
	except ValueError:
		sys.stderr.write( "ValueError: fatal error running blast.\n" )
		sys.exit(1)
	except IOError:
		sys.stderr.write( "IOError: fatal error running blast.\n" )
		sys.exit(1)
	if retcode == 13:
		sys.stderr.write( "Permission Denied Error: fatal error running blast." 
          "Is the blastall binary in the path with execution and read permissions?\n" )
		sys.exit(1)			
## End of exe_blast

def extractPE( infile, outfile, seqfile, alignLength, alignIdentity ):
	ifh=open(infile,'r')
	ofh=open(outfile,'w')
	tags={}
	while 1:
		line=ifh.readline().rstrip('\n')
		if not line:
			break
		col=line.split('\t')
		if (float(col[2]) < alignIdentity or int(col[3]) < alignLength):
			continue
		tag=re.search('(.+)\/\d{1}$',col[0]).group(1)
		if tag not in tags:
			tags[tag]=1
		tags[tag]+=1
	ifh.close()
	
	sfh=open(seqfile,'r')
	while 1:
		tagA=sfh.readline().rstrip('\n')
		if not tagA:
			break
		seqA=sfh.readline().rstrip('\n')
		tagB=sfh.readline().rstrip('\n')
		seqB=sfh.readline().rstrip('\n')
		tag=re.search('\>(.+)\/\d{1}$',tagA).group(1)
		if tag in tags and tags[tag]==2:
			ofh.write(tagA+'\n'+seqA+'\n'+tagB+'\n'+seqB+'\n')
	sfh.close()
	ofh.close()
## End of extractPE

def findOTU(desc, taxLevel):
	col=desc.split(';')
	for c in col:
		ele=c.split('__')
		if ele[0] == taxLevel:
			return ele[1]
## End of findOTU
	
def taxoAssign(infile, outfile, libfile, taxLevel, alignLength, alignIdentity):
	ifh=open(infile,'r')
	ofh=open(outfile,'w')
	lfh=open(libfile,'r')
	
	lib={}
	while 1:
		line=lfh.readline().strip('\n')
		if not line:
			break
		col=line.split('\t')
		lib[col[0]]=col[1]
	lfh.close()
	
	taxRes={}
	while 1:
		line=ifh.readline().strip('\n')
		if not line:
			break
		col=line.split('\t')
		if (float(col[2]) < alignIdentity or int(col[3]) < alignLength):
			continue
		tag, end = re.search('(.+)\/(\d{1})$',col[0]).group(1,2)
		desc=lib[col[1]]
		otu=findOTU(desc,taxLevel)
		
		if tag not in taxRes:
			taxRes[tag] = [None, None]
		
		if end == '1':
			taxRes[tag][0] = otu
		elif end == '2':
			taxRes[tag][1] = otu
	ifh.close()
		
	for tag in taxRes:
		if (taxRes[tag][0] and taxRes[tag][1]) and (taxRes[tag][0] == taxRes[tag][1]):
			ofh.write(tag+'\t'+taxRes[tag][0]+'\n')
	ofh.close()
## End of taxoAssign

def taxoStat(infile, seqfile, statfile, readfile):
	ifh=open(infile,'r')
	sfh=open(seqfile,'r')
	statfh=open(statfile,'w')
	rfh=open(readfile,'w')
	
	otus={}
	tags={}
	while 1:
		line=ifh.readline().rstrip('\n')
		if not line:
			break
		col=line.split('\t')
		if col[1] not in otus:
			otus[col[1]] = 0
		otus[col[1]]+=2
		tags[col[0]]=col[1]
	ifh.close()
	
	while 1:
		tagA=sfh.readline().rstrip('\n')
		if not tagA:
			break
		seqA=sfh.readline().rstrip('\n')
		tagB=sfh.readline().rstrip('\n')
		seqB=sfh.readline().rstrip('\n')
		
		tag=re.search('\>(.+)\/\d{1}$',tagA).group(1)
		if tag in tags:
			otu=tags[tag]
			newTagA=tagA+'|'+otu
			newTagB=tagB+'|'+otu
			rfh.write(newTagA+'\n'+seqA+'\n'+newTagB+'\n'+seqB+'\n')
	sfh.close()
	rfh.close()
	
	for otu in otus:
		statfh.write(otu+'\t'+str(otus[otu])+'\n')
	statfh.close()
			
#=================================#

if __name__=='__main__':
	pars = read_params( sys.argv )

	if pars['version']:
		sys.stdout.write("MG16S version "+__version__+"\t("+__date__+")"+"\n")
		sys.exit(0)
		
	if pars['in'] is None: 
		sys.stderr.write( "The INPUT_FILE parameter need top be specified\n"
                          "Type MG16S.py -h for more info\n")
		sys.exit(0)
		
	if pars['outprefix'] is None:
		pars['outprefix']='MG16S_output'
		sys.stdout.write( "[Warning]: You did not provide a prefix for output files,\n"
						"'MG16S_output' is used as default output prefix\n")
	if pars['probe'] is None:
		sys.stderr.write( "The --probe parameter need top be specified\n"
                          "Type MG16S.py -h for more info\n")
		sys.exit(0)
		
	if pars['GreenGenesDB'] is None:
		sys.stderr.write( "The --GreenGenesDB parameter need top be specified\n"
                          "Type MG16S.py -h for more info\n")
		sys.exit(0)
		
	if pars['GreenGenesTax'] is None:
		sys.stderr.write( "The --GreenGenesTax parameter need top be specified\n"
                          "Type MG16S.py -h for more info\n")
		sys.exit(0)
		
	## Print out selected option
	printOptions(pars)
		
	## start run infile vs. probe blastn search
	sys.stdout.write("## Now fishing potential 16S reads using the probe sequence...\n")
	run_blastn(pars['formatdb'],pars['blast'],pars['in'],pars['outprefix']+'.Probe.bls',
				pars['probe'],pars['evalue'],pars['nproc'])
				
	## extract PE reads with both ends matching probe 16S (extractPaired.pl and extractSeq.pl)
	sys.stdout.write("## Extracting PE reads with both ends matched to the probe sequence...\n")
	extractPE( pars['outprefix']+'.Probe.bls', pars['outprefix']+'.PE.fa', 
				pars['in'], int(pars['alignLength']), 70)
	
	## run extracted reads against greengenesDB
	sys.stdout.write("## Extraction done!\n"
					"## Now running BLASTN search against GreenGenes 16S database...\n")
	run_blastn(pars['formatdb'], pars['blast'], pars['outprefix']+'.PE.fa',
				pars['outprefix']+'.GreenGenes.bls', pars['GreenGenesDB'],
				pars['evalue'], pars['nproc'])
	
	## assign taxonomy information to blastn output (taxoAssign.pl)
	sys.stdout.write("## Taxonomy information integration...\n")
	taxoAssign( pars['outprefix']+'.GreenGenes.bls', pars['outprefix']+'.GreenGenes.Anno.txt',
				pars['GreenGenesTax'], pars['taxLevel'], 
				int(pars['alignLength']), float(pars['alignIdentity']))
	
	taxoStatOutfile=pars['outprefix']+'.'+tax_lookup[pars['taxLevel']]+'.stat.txt'
	readOutfile=pars['outprefix']+'.'+tax_lookup[pars['taxLevel']]+'.reads.fa'
	taxoStat( pars['outprefix']+'.GreenGenes.Anno.txt', pars['in'], taxoStatOutfile, readOutfile)
	
	sys.stdout.write("## Task finished, bye!\n")	
	exit(0);
	
## END
