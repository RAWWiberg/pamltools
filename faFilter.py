#!/usr/bin/env python3
"""
This script takes a fasta file and removes samples that occur in a list
requires biopython
"""

import sys
import argparse
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord

def fasfilter(fas,to_drop):
    fastafile = SeqIO.parse(open(fas,'r'),"fasta")
    for sequence in fastafile:
        seqid = sequence.id
        if seqid in to_drop:
            continue
        else:
            print(sequence.format("fasta-2line"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-fas',
                       help = 'fasta format file for checking')

    parser.add_argument('-t',
                       help = 'list of species to drop')
    
    args = vars(parser.parse_args())

    try:
        fas = args['fas']
        dropfil = open(args['t'], 'r')
        to_drop = dropfil.readlines()
        to_drop = [i.replace('\n','') for i in to_drop]
        #print(to_drop)
        fasfilter(fas,to_drop)
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass


