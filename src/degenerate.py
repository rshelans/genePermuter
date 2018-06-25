#! /usr/bin/env python2.7
import sys
import collections
import argparse
import re
import itertools
from Bio import Seq
import random as rd	

def main(args):
	k=6
	BASE=set(['A','C','G','T'])
	perms=["".join(z) for z in itertools.product(BASE,repeat=k)]
	for m in perms:
		if m.count('TGA')>0 or m.count('TAA')>0 or m.count('TAG')>0:	
			print "{}\t{}".format(m,str(Seq.Seq(m).reverse_complement()))
	#keep=True
	#for perm in perms:
	  #keep = True
	  #if perm.count('Y') ==1:
	  	#d=degenerate(perm)
	  	#for m in d:
	  		#if m.count('TGA')==0 and m.count('TAA')==0 and m.count('TAG')==0:	
	  			#pass
	  		#else:
	  			#keep=False
	  	#if keep:		
	  		#print "{}\t{}".format(perm,"\t".join([x for x in d]))
		
	#for x in itertools.permutations(list(args[1])):
		#kmer="".join(x)
		#print "{}\t{}".format(kmer,"\t".join(degenerate(kmer)))


	#print "c({})".format(",".join(["'"+x+"'" for x in degenerate(args[1])]))
	#for m in perms :
	#	k=[]
	#	for d in degenerate(m):
	#		k.append(d)
		  #k.append(str(Seq.Seq(d).reverse_complement()))
	#	print "{}\t{}".format(m,"\t".join(k))
	#print "c({})".format(",".join(["'"+x+"'" for x in k]))




def degenerate(word):
   d=[word] #of degenerates formed from word 
   deg={'W':['A','T'],'S':['C','G'],'M':['A','C'],'K':['G','T']\
       ,'R':['A','G'],'Y':['C','T'],'B':['C','G','T']\
       ,'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G']\
       ,'N':['A','T','G','C']}
   for i in range(len(word)):
       if word[i] in deg:
           n=list()
           for word in d:
               n+=[word.replace(word[i],x,1) for x in deg[word[i]]]
           d=n
   return d

if __name__ == "__main__" :
    try:
        sys.exit(main(sys.argv))
    except EnvironmentError as (errno,strerr):
        sys.stderr.write("ERROR: " + strerr + "\n")
        sys.exit(errno)
        
        
        
