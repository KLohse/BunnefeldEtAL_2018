#!/usr/bin/env python
from __future__ import division 
import itertools, sys
import numpy as np

def piMaker(blockCounts,samps,block_len):
	Test=[[(iii[samps[0]],iii[samps[1]]) for iii in block] for block in blockCounts]
	Test2=[[1 for iii in block if iii[0] != iii[1]] for block in Test]
	Test3=[sum(block)/block_len for block in Test2]
	return(float(sum(Test3))/len(Test3))

def main():
	myBlockDict=np.load(sys.argv[1]).item()
	blockCounts = myBlockDict.values()#Ignoring
	block_len=int(sys.argv[2])
	sam=[int(i) for i in list(sys.argv[3])]
	samPairs=[]
	for subset in itertools.combinations(sam, 2):
		samPairs.append(list(subset))

	pi=[]
	for i in range(len(samPairs)):
		pi1=piMaker(blockCounts, samPairs[i],block_len)
		pi.append(pi1)	
		#print(pi1)
	print(pi)
if __name__== "__main__":
	main()
