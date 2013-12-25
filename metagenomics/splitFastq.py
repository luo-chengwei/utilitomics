import sys
import os

ifh1 = open(sys.argv[1],'r')
ifh2 = open(sys.argv[2],'r')
chunk = int(sys.argv[3])
prefix = sys.argv[4]

ofhs = []
for i in range(chunk):
	ofh1 = open(prefix+'.1.'+str(i+1)+'.fastq', 'w')
	ofh2 = open(prefix+'.2.'+str(i+1)+'.fastq', 'w')
	ofhs.append((ofh1, ofh2))

index=0
while 1:
	tag1=ifh1.readline().rstrip('\n')
	if not tag1:
		break
	seq1=ifh1.readline().rstrip('\n')
	dis1=ifh1.readline().rstrip('\n')
	qual1=ifh1.readline().rstrip('\n')
	tag2=ifh2.readline().rstrip('\n')
	seq2=ifh2.readline().rstrip('\n')
	dis2=ifh2.readline().rstrip('\n')
	qual2=ifh2.readline().rstrip('\n')
	index+=1
	res = index%chunk
	ofhs[res][0].write('%s\n%s\n%s\n%s\n' % (tag1, seq1, dis1, qual1))
	ofhs[res][1].write('%s\n%s\n%s\n%s\n' % (tag2, seq2, dis2, qual2))
ifh1.close()
ifh2.close()

for ofh1, ofh2 in ofhs:
	ofh1.close()
	ofh2.close()

