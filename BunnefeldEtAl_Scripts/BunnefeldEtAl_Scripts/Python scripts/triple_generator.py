#!/usr/bin/env python
import sys, itertools, vcf
import numpy as np

#this creates a dictionary per block, sums the mutation types per block and spits out the counts as a list  
def countTypes(l,mutTypes):
	myCounter={}
	for i in range(0,len(mutTypes)):
     		myCounter.update({mutTypes[i]:0})
	if len(l) > 0:
		for j in range(len(l)):
	        	myCounter[l[j]] += 1
	return(RootlessDict(myCounter).values())

#same function but exports the dict (so you can check the order of the counts) I am not sure if these are fixed...need to check for another species
def countTypesD(l,mutTypes):
	myCounter={}
	for i in range(0,len(mutTypes)):
     		myCounter.update({mutTypes[i]:0})
	if len(l) > 0:
		for j in range(len(l)):
	        	myCounter[l[j]] += 1
	return(RootlessDict(myCounter))
##this is the order...{('0', '0', '1'): 5, ('1', '0', '1'): 0, ('1', '1', '1'): 0, ('1', '0', '0'): 1}

#unrooting the mutation type counts (without a root config 001 is the same 110, etc.)
def RootlessDict(MutCountsDict):		    
	myNewCounter={}
	for it in itertools.combinations(MutCountsDict.keys(),2):
		if ([int(a)+int(b) for a,b in zip(it[0],it[1])]==[1 for x in range(len(MutCountsDict.keys()[0]))]):
			myNewCounter.update({it[0]:MutCountsDict[it[0]]+MutCountsDict[it[1]]})
	return(myNewCounter)

#This is a function to extract counts of mutation types per triple. It requires a list of dictionary values where each list is a block with polymorphic sites represented as tuples.
#It also requires a three element list indicating which samples to take (usually from a choice of six)
#muttypes are the possible mutation types for a triple sample and myorder is reordering vector needed after you have unrooted the block (this is hardcoded right now, might want to work out how to generalise this)

def tripleMaker(blockCounts,samps,mutTypes,myorder):
	yo3=[[(iii[samps[0]],iii[samps[1]],iii[samps[2]]) for iii in block] for block in blockCounts]
	countYo=[countTypes(pp,mutTypes) for pp in yo3]
	countYo2 = [[mylist[i] for i in myorder] for mylist in countYo]
	return(countYo2)	
	
def main():
	myBlockDict=np.load(sys.argv[1]).item()#filtered set of blocks
	blockCounts = myBlockDict.values()
	tripleSet=[[0,1,2],[0,1,3],[0,1,4],[0,1,5],[2,3,4],[2,3,5],[0,2,3],[1,2,3],[0,4,5],[1,4,5],[2,4,5],[3,4,5],[0,2,4],[1,3,5]]#this could be a sys.argv, esp. for species with less than six samples
	myorder=[3,1,0,2]#new order is 100,010,001,111
	mutTypes=[it for it in itertools.product('01',repeat=3)]
	AllTriples=[tripleMaker(blockCounts,tripleSet[i],mutTypes,myorder) for i in range(len(tripleSet))]
	BlocksAllT=[]
	BlocksAllT.append(myBlockDict.keys())
	BlocksAllT.append(AllTriples)
	out=str(BlocksAllT).replace("[","{").replace("]","}").replace(")","}").replace("(","{").replace("'","\"")
	print(out)
if __name__== "__main__":
	main()
