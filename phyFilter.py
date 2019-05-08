#!/usr/bin/env python3
"""
This script takes a phylip file and removes samples that occur in a list

"""

import sys
import argparse

def phyfilter(phy,to_drop):
    alignments = open(phy, 'r')
    content = alignments.readlines()
    # reads in contents of phylip file to a list, if file is too large
    # this might be a problem
    content = [i.replace('\n', '') for i in content]
    head = content[0].split()
    #  
    i=1
    while i < len(content):
        linedat=content[i].split()
        ID=linedat[0]
        if ID in to_drop:
            del content[i]
            head=[str(int(head[0])-1),head[1]]
            content[0]="    ".join(head)
        else:
            i=i+1
    for i in range(0,len(content)):
        print(content[i])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-phy',
                       help = 'phylip format file for checking')

    parser.add_argument('-t',
                       help = 'list of species to drop')
    
    args = vars(parser.parse_args())

    try:
        phy = args['phy']
        dropfil = open(args['t'], 'r')
        to_drop = dropfil.readlines()
        to_drop = [i.replace('\n','') for i in to_drop]
        #print(to_drop)
        phyfilter(phy,to_drop)
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass


