#!/usr/bin/env python
from __future__ import division
import math, sys
import numpy as np

#this outputs a table of counts of singletons, doubletons etc for each block

#./nton_counts_table.py inputfile.npy number_of_indivs > Cfun_nton_counts_table_output.txt &

def nton_counts(d,n):
	newlist = []
	for bl in d.values():
		emptlist = [0]*(n-1)
		for site in bl:
			tonnage = site.count('1')	#counts '1's in each site
			emptlist[(tonnage-1)] += 1	#adds the count to emptlist
		blntc = [(emptlist[i] + emptlist[-i-1]) for i in range(int(math.ceil(len(emptlist) / 2)))]	#folds counts by adding first and last, second and second last etc.
		if len(emptlist)%2 == 1:
			blntc[-1] = int(blntc[-1] / 2)	#divides last count by 2 if the number of individuals is even (so the middle one is not counted twice)
		newlist.append(blntc)
	newdict = dict(zip(d.keys(),newlist)) #turns it into a dictionary with blocks as keys and counts as values
#	header = ['block','singletons','doubletons','tripletons','quadrupletons','quintupletons'][0:(int(math.floor(n/2))+1)] #header of table, cut to length required (works for up to ten individuals)
	return newdict

def main():
	myBlockDict = np.load(sys.argv[1]).item()
	countdict =  nton_counts(myBlockDict,int(sys.argv[2]))
	outlist = []
	for key, value in countdict.iteritems():
		outlist.append([key, value])
	out=str(outlist).replace("[","{").replace("]","}").replace(")","}").replace("(","{")
	print out
if __name__== "__main__":
	main()
