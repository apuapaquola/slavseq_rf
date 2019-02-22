#!/usr/bin/env python
'''This program takes headerless paired-end SAM input sorted by read name. For a read pair, the main alignment is defined as the primary alignment of R1 (i.e., the alignment not having flags 2048 or 256). The coordinates of that alignment are used to index all other alignments belonging to the pair (i.e. R2 and secondary alignments) using a tabix index.

Example:
id1/r1: maps to chr1:1000-1080 (primary r1 alignment)
id1/r1: maps to chr14:2020-2040 (secondary r1 alignment)
id1/r2: maps to chrX:10-110 (primary r2 alignment)

This way, a search for chr1:1000-2000 will return all 3 alignments. A
search for chr14:2020-2040 or chrX:10-110 will return nothing.
'''

import sys
import re

def is_primary_r1_alignment(a):
    return int(a[1]) & (2048|256|4) == 0 and int(a[1]) & (64) != 0
    
def qname(a):
    m=re.search('^([^/]*)', a[0])
    return m.group(1)

def cigar2reflen(cigar):
    r=0
    for m in re.finditer('(\d+)[MDN]',cigar):
        r+=int(m.group(1))
    return r

def output_alignments(alignments, primary_r1):
    if primary_r1 is None and len(alignments)>0:
        print(alignments,file=sys.stderr)
    if len(alignments)>0 and primary_r1 is not None:
        for al in alignments:
            print('\t'.join([primary_r1[2], primary_r1[3], str(int(primary_r1[3])+cigar2reflen(primary_r1[5])), al]),end='')
            
            #print(alignments)


alignments=[]
primary_r1=None

last_qname=None

for line in sys.stdin:
    a=line.split('\t')
    a[-1]=a[-1].rstrip()

    current_qname=qname(a)
    
    if last_qname != current_qname:
        output_alignments(alignments, primary_r1)
        alignments=[]
        primary_r1=None
        last_qname=current_qname

    alignments.append(line)
    #alignments.append(str(cigar2reflen(a[5]))+'_'+str(((int(a[1])&128)//128)+1))

    if is_primary_r1_alignment(a):
        if primary_r1 is not None:
            print("Primary R1 error:", primary_r1[0], file=sys.stderr, flush=True)
            assert primary_r1 is None
            
        primary_r1=a

output_alignments(alignments, primary_r1)
