#usage python drawMapping.py [gene fasta] [bls file] [out figure]
import sys
import matplotlib.pyplot as plt

ifh=open(sys.argv[1],'r')
bfh=open(sys.argv[2],'r')
outfig=sys.argv[3]
xlabel='Position on 16S rRNA (bp)'
ylabel='Identity'

region_min=[69,137,433,576,822,986,1117,1243,1435]
region_max=[99,242,497,682,879,1043,1173,1294,1465]


y=[]
for i in region_min:
	y.append(100)

lib=''
while 1:
	line=ifh.readline().rstrip('\n')
	if not line:
		break
	if line.count('>')!=0:
		continue
	lib+=line
ifh.close()
length=len(lib)

id=[]
xmin=[]
xmax=[]
while 1:
	line=bfh.readline().rstrip('\n')
	if not line:
		break
	col=line.split('\t')
	start=float(col[8])
	end=float(col[9])
	if start>end:
		start,end=end,start
	identity=float(col[2])
	id.append(identity)
	xmin.append(start)
	xmax.append(end)
bfh.close()

fig=plt.figure(figsize=(32,8))
ax=fig.add_subplot(111)
ax.hlines(id, xmin, xmax, colors=(0,0,0,.5), lw=1, alpha=0.2)
ax.hlines(y, region_min, region_max, colors='r', lw=5)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(0,length)
ax.set_ylim(60,100)
fig.savefig(outfig,format='eps')