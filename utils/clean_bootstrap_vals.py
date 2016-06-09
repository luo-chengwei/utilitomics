import re, os, sys

tree = open(sys.argv[1], 'r').read()
new_tree = re.sub('\)\d{1}\.\d{3}', ')', tree)
ofh = open(sys.argv[2], 'w')
ofh.write(new_tree+'\n')
ofh.close()
