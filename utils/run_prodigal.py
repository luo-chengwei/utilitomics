import sys, re, os
from Bio import SeqIO
from subprocess import PIPE, call, Popen, STDOUT
import multiprocessing as mp

DEVNULL = open(os.devnull, 'wb')

def run_prodigal(args):
	prefix, mode, infile = args
	call(['prodigal', '-a', prefix+'.faa', '-d', prefix+'.ffn',
			 '-i', infile, '-o', prefix+'.gbk', '-p', mode], 
			 stdout = DEVNULL, stderr = STDOUT)
	return 0

def chunks(l, n):
	for i in range(0,len(l), n):
		yield l[i:i+n]

def main():
	try:
		infile, outprefix, mode, n_st = sys.argv[1:]
		ncpu = int(n_st)
	except:
		print "Usage: python run_prodigal.py [infile] [out prefix] [meta|single] [ncpu]"
		exit(1)
	# split infile
	p1 = Popen(['grep', '>', infile], stdout = PIPE, stderr = sys.stderr)
	p2 = Popen(['wc', '-l'], stdin = p1.stdout, stdout = PIPE, stderr = sys.stderr)
	o, e = p2.communicate()
	numRec = int(o.replace(' ', '', 100))
	n = max(1, numRec/ncpu)
	X = list(chunks(range(numRec), n))
	if len(X) > ncpu:
		A = [item for l in X[ncpu-1:] for item in l]
		X = X[:ncpu-1]+[A]
	temp_infiles = []
	for i in range(ncpu): temp_infiles.append('%s.input_%i' % (outprefix, i))
	ofhs = []
	for temp_infile in temp_infiles: ofhs.append(open(temp_infile, 'wb'))
	Y = {}
	for ofh_ind, inds in enumerate(X):
		for ind in inds: Y[ind] = ofh_ind
	for rec_ind, rec in enumerate(SeqIO.parse(infile, 'fasta')):
		ofh_ind = Y[rec_ind]
		ofhs[ofh_ind].write('>%s\n%s\n' % (rec.description, rec.seq))
	for ofh in ofhs: ofh.close()
	
	# make mp cmds
	cmds = []
	for temp_infile in temp_infiles:
		sub_outprefix = temp_infile
		cmds.append([sub_outprefix, 'single', temp_infile])
	
	pool = mp.Pool(ncpu)
	pool.map_async(run_prodigal, cmds)
	pool.close()
	pool.join()
	
#	for cmd in cmds: run_prodigal(cmd)

	# merge_results
	faa_fh = open('%s.faa' % outprefix, 'wb')
	ffn_fh = open('%s.ffn' % outprefix, 'wb')
	for i in range(ncpu):
		temp_faa = '%s.input_%i.faa' % (outprefix, i)
		temp_ffn = '%s.input_%i.ffn' % (outprefix, i)
		temp_gbk = '%s.input_%i.gbk' % (outprefix, i)
		temp_infile = '%s.input_%i' % (outprefix, i)
		# convert faa
		for rec in SeqIO.parse(temp_faa, 'fasta'):
			desc = rec.description.split(' # ')
			tag = re.sub('\_\d+$', '', desc[0])
			start, end, strand = desc[1:4]
			if strand == '1': strand = '+'
			else: strand = '-'
			geneID = desc[4].split(';')[0]
			faa_fh.write('>%s|%s|%s-%s|%s\n%s\n' % 
				(tag, geneID, start, end, strand, 
				str(rec.seq).replace('*', '', 100)))
		
		# convert ffn	
		for rec in SeqIO.parse(temp_ffn, 'fasta'):
			desc = rec.description.split(' # ')
			tag = re.sub('\_\d+$', '', desc[0])
			start, end, strand = desc[1:4]
			if strand == '1': strand = '+'
			else: strand = '-'
			geneID = desc[4].split(';')[0]
			ffn_fh.write('>%s|%s|%s-%s|%s\n%s\n' % 
				(tag, geneID, start, end, strand, 
				str(rec.seq).replace('*', '', 100)))
		os.unlink(temp_gbk)
		os.unlink(temp_faa)
		os.unlink(temp_ffn)
		os.unlink(temp_infile)
		
		
	faa_fh.close()
	ffn_fh.close()		
	
if __name__ == '__main__': main()
