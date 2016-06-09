import sys
import os
import re
from subprocess import PIPE, Popen, call

fq1, fq2, db, prefix = sys.argv[1:]
bowtie2_logfh = open(prefix+'.bowtie2.log','w')
bamfile = prefix+'.bam'
bowtie2_cmd = ['bowtie2', '-x', db, '-1', fq1, '-2', fq2]
samtools_view = ['samtools', 'view', '-bhS', '-']
samtools_sort = ['samtools', 'sort', '-', prefix]
samtools_index = ['samtools', 'index', bamfile]
p1 = Popen(bowtie2_cmd, stdout = PIPE, stderr = bowtie2_logfh)
p2 = Popen(samtools_view, stdin = p1.stdout, stdout = PIPE, stderr = bowtie2_logfh)
p3 = Popen(samtools_sort, stdin = p2.stdout, stdout = PIPE, stderr = bowtie2_logfh)
p1.stdout.close()
p2.stdout.close()
output, err = p3.communicate()
samtools_index = ['samtools', 'index', bamfile]
call(samtools_index, stderr = bowtie2_logfh, stdout = bowtie2_logfh)
bowtie2_logfh.close()


