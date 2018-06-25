#! /usr/bin/env python2.7
'''
2013
@author: Robert Shelansky
Description:

'''
import sys
import argparse
from Bio import SeqIO
import collections
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
def main(args):
    """
    Does All the grunt work
    """
    opts = parse_args()

    remove=[]
    if opts.remove is not None: 
        for line in open(opts.remove,'r'):
            remove.append(line.strip())
    fasta_id  = {}
    name_id   = collections.defaultdict(lambda:'')
    for i,record in enumerate(SeqIO.parse(sys.stdin,"fasta")):
       fasta_id[record.id]=record
       name_id[record.description.split()[1].strip()]=record.id
    total=len(fasta_id.keys())
    not_removed=[]
    #c1=0
    #c2=0
    for r in remove:
        if fasta_id.get(name_id[r]) is not None:       
            del fasta_id[name_id[r]]
            #c1+=1
        elif fasta_id.get(r) is not None:
            #c2+=1
            del fasta_id[r]
        else:
            not_removed.append(r)
    removed=len(fasta_id.keys())
    SeqIO.write([fasta_id[i] for i in fasta_id.keys()],sys.stdout,"fasta")
    sys.stderr.write("{}\n".format(opts.remove))
    sys.stderr.write("Removed {} from {} expected {} Missed {}.\n".format(total-removed,total,len(remove),len(not_removed)))

    sys.stderr.write("\n".join(not_removed)+"\n")
    #sys.stderr.write("{},{}\n".format(c1,c2))
def parse_args():
    parser = argparse.ArgumentParser(description=\
"""
""")

    parser.add_argument("-r","--remove",type=str,default = None)
    args = parser.parse_args()
    return args
#MAKE COMMAND LINE RUNNABLE
if __name__ == "__main__" :
    #try:
        sys.exit(main(sys.argv))
    #except EnvironmentError as (errno,strerr):
    #    sys.stderr.write("ERROR: " + strerr + "\n")
    #    sys.exit(errno)
#****************************************************************************

