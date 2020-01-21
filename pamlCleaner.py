#!/usr/bin/env python3
"""
This script takes a phylip file and checks that the formatting is
appropriate for paml:

Checks sequence name length.
Checks length divisible by 3 (i.e. there is a full ORF)

if all above checks pass:
Also removes terminal STOP codons. 
"""

import sys
import argparse
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Gapped, generic_dna
import re
import os


def pamlcleaner(phy,min_aln_l):
    output=phy.replace('.phy', '_clean.phy')
    outputfile = open(output, 'w')
    alignments = open(phy, 'r')
    content = alignments.readlines()
    # reads in contents of phylip file to a list, if file is too large
    # this might be a problem
    content = [i.replace('\n', '') for i in content]
    head = content[0].split()
    # Check if alignment length is > min_aln_l
    if int(head[1]) > min_aln_l:
        seqs = [] # everything in file except header
        ids = []
        for line in content[1:]:
            if line != '':
                ID=line.split()[0]
                s=line.split()[1]
                ids.append(ID)
                seqs.append(s)
        # convert each string into a seq object which can be checked by
        # BioPython functions
        seqs = [SeqRecord(Seq(seq, alphabet = Gapped(generic_dna, "-"))) for seq in seqs]
        if len(seqs) != int(head[1]):
            # number of sequences does not match header
            sys.stderr.write("more sequences than mentioned in file header: "+
                             phy+" EXITING!\n")
            # Remove the file and exit
            os.remove(output)
            sys.exit()
        name = 0
        end_stop=0
        for s in seqs:
            if str(s.translate(1))[:-1].count("*") > 0:
                # premature stop codons found
                sys.stderr.write("premature stop codon found in: "+phy+
                                 " - "+ids[name]+" EXITING!\n")
                                    # remove the file and exit
                os.remove(output)
                sys.exit()
            if len(s.seq) != int(head[2]) or (int(head[2])-3)%3 != 0:
                # sequence not same length as described in header
                sys.stderr.write("sequence lengths do not correspond"+
                                 " to header OR not a multiple of 3 "+
                                 ": "+phy +" - "+ids[name]+" EXITING!\n")
            if str(s.translate(1))[:-1][-1] == '*':
                    # last codon is stop codon
                    end_stop=1
            name=name+1
        # all checks have passed, print the sequences (minus
        # last three bases (terminal STOP codon))
        if end_stop == 1:
            # Write the sequence minus the last three bases for each codon.
            head[2] = str(int(head[2])-3)
            outputfile.write(' '.join(head)+'\n')
            for i in range(0,len(ids)):
                outputfile.write(ids[i]+'  '+str(seqs[i].seq)[:-3]+'\n')
        else:
            outputfile.write(' '.join(head)+'\n')
            for i in range(0,len(ids)):
                outputfile.write(ids[i]+'  '+str(seqs[i].seq)+'\n')
    else:
        sys.stderr.write("Alignment shorter than : "+str(min_aln_l)+" - "+
                             phy+" EXITING!\n")
        os.remove(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-phy',
                       help = 'phylip format file for checking')

    parser.add_argument('-min_aln_l',default=100,type=int,
                       help = 'minimum alignment length')
    
    args = vars(parser.parse_args())

    try:
        phy = args['phy']
        pamlcleaner(phy,args['min_aln_l'])
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass


