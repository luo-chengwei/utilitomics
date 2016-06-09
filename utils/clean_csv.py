import sys, re, os

ifh = open(sys.argv[1],'r')
ofh = open(sys.argv[1]+'.tmp','w')
file = ifh.read()
lines = file.split(chr(13))
for line in lines: ofh.write(line+'\n')
ifh.close()
ofh.close()
os.rename(sys.argv[1]+'.tmp', sys.argv[1])
