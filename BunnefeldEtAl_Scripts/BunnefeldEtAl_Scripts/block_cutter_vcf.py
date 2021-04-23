#!/usr/bin/env python
import sys, itertools
import vcf
import numpy as np

#./block_cutter_vcf_LB.py Cfun.CallableLoci.mbq10mmq20md1.masked.CALLABLE.bed2 Cfun.calibrated_haploid.filtered_variants.CallableLoci.md1.masked.vcf.gz loci.txt 271 Cfun135.Cfun142.Cfun83.Cfun71.Cfun86B.Cfun3549B 
"""
Updates a dictionary by appending tuples
"""
def update_dict_list(d, l):
	ke = str(l[0])
	vv = map(int,l[1:])
	if d.has_key(ke): 
		d[ke].append(vv)
	else:
		d.update({ke:[vv]})
	return d

"""
Parses a bed file line by line into a dict (using update_dict_list)
"""
def bed_to_dict(bed):
	ldict={}
	with open(bed) as f:
		for line in f:
			l = line.rstrip("\n").split("\t")
			ll = update_dict_list(ldict, l)
	return ll

"""
Filters dict by the maximum length post-filtering
"""
def filter_dict_maxcov(d, le):
	dnew = {k:v for k,v in d.iteritems() if sum([item[-1] for item in map(list,v)]) > le}
	return dnew

"""
Filters dict by the maximum length each block spans pre-filtering.
"""
def filter_dict_maxspan(d,le):
	dnew = {k:v for k,v in d.iteritems() if v[-1][1]-v[0][0] <= le}
	return dnew

"""
Adds coordinates post filtering. These are used by blockcut.
"""
def	fltrd_coord(ttl): 
	count= 0
	nl=[]
	for i in ttl:
		count = count + i[-1]
		nl.append(i[0:2]+[count])
	return nl

"""
Takes an ordered list of range coordinates for a contig and returns a nested list of range coordinates for blocks of length bl within that contig. Ranges in the original list that span twor or more blocks are split up in the process.   
"""
def blockcut(inlist, bl):
	fin_block_list = []
	left_list = inlist
	max_len = inlist[-1][-1]
	for i in range(bl,max_len,bl):
		new_bl_list = list(itertools.takewhile(lambda v: v[-1] <= i, left_list))
		left_list = list(itertools.dropwhile(lambda v: v[-1] <= i, left_list))
# this deals with the rare case of ranges that have exactly the length of a block:
		if len(new_bl_list) == 1 and new_bl_list[0][-1] == i:
			pass
# a block boundary spans a range, this is the most likely case:
		else:
			nb = left_list[0]
#we need to break the range of the first block in left_list to make a new coordinate which is used to make two new blocks: new_last block and new fisrt block.
			if len(new_bl_list)>0:
				new_coor = nb[1]-nb[-1]+i
			else:			
				new_coor = nb[0]+bl
#recut the last new block and append it to new_bl_list, append new_bl_list to fin_block_list
			new_last_block = [nb[0],new_coor,i]
			new_bl_list.append(new_last_block)			
#recut the first left block and append it to left_list, append new_bl_list to fin_block_list
			new_first_block = [new_coor,nb[1],nb[-1]]
			left_list.pop(0)
			left_list[:0] = [new_first_block]
		fin_block_list.append(new_bl_list)
	return fin_block_list

"""
Applies blockcut to a dict of chrom:coordinates and makes a new dict of (chrom, blocksstart, blockend): coordinates
"""
def blockcutdict(tdict, bl):
	newdict={}
	for ke, vv in tdict.iteritems():
		vvnew =  blockcut(fltrd_coord(vv), bl)
		for v in vvnew:
			knew = tuple([ke,v[0][0],v[-1][1]])
			newdict.update({knew:v})
	return newdict

"""
Takes a column from a vcf file and extract variable sites as a tuple. At the moment GQ and DP tresholds are not used; this could be added easily (see ness.vcf)
"""
def get_vcf_column(record, min_GQ=0, min_DP=0, samples='all'):
	GTs = []
	if record.is_snp:
		if samples == 'all' or samples == ['all']:
			for sample in record.samples:
				GTs.append(sample['GT'])
		else:
			for i in samples:
				GTs.append(record.genotype(i)['GT'])
	return tuple(GTs)

"""
Opens a region from a vcf file and extracts variable sites as a list of tuples. At this stage only sites that are biallelic and contain genotype info for all samples are retained.
"""
def get_vcf_region(vcf_reader, chrom, region, min_GQ=0, min_DP=0, samples='all'):
	# Expects samtools style region and a tabixed vcf "chromosome:start-end"
	# Note that i) the tabixed vcf has to end in .gz rather than .bgz for some buggy reason and ii) bgzip needs to be used for the zipping.
	start_coord = region[0]
	end_coord = region[1]
	ll = []
	for record in vcf_reader.fetch(chrom,start_coord,end_coord):
		if record.is_snp:
			ttu = get_vcf_column(record, min_GQ, min_DP, samples)
	#Filters out sites that are not variable in the ingroup:
			if ttu.count(ttu[0]) == len(ttu):
	#or None in ttu:
				pass
			else:
				ll.append(ttu)
	return ll

def get_vcf_contigs(vcf_reader, tdict, chromset, min_GQ=0, min_DP=0, samples='all'):
	newdict={}
	for ke, vv in tdict.iteritems():
		chrom = ke[0]
		blockl = []
		if ke[0] in chromset:
			for i in vv:
				blockl.extend(get_vcf_region(vcf_reader, chrom, i[0:2], min_GQ, min_DP, samples))
			newdict.update({ke:blockl})
		else:
			for i in vv:
				blockl.extend([])
			newdict.update({ke:blockl})
	return newdict	

def get_chrom_set(fil):
	chromset = set()
	with open(fil) as f:
		for line in f:
			chromset.add(line.rstrip("\n"))
	return chromset

###THIS REMOVES BAD BLOCKS (i.e. ones with more than le Nones (e.g. 5)
def filter_dict_ByNoneCount(d,le):
	dnew = {k: v for k, v in d.iteritems() if sum(iii.count(None) for iii in v) <= le}
	return dnew

###THIS IS LESS EXTREME, BUT REMOVES SITES FROM BLOCKS WHICH HAVE NONES OR MULTIALLELIC SITES (i.e., makes them invariant)
def change_dict_RemoveMultiSNoneS(d):
	for k,v in d.iteritems():
		v2=[iii for iii in v if (iii.count('2') == 0) and (iii.count('3') == 0) and (iii.count('4') == 0) and (iii.count(None) == 0)]
		d.update({k:v2})
	return d

"""
This function summarises blockwise data (a dictionary) from two populations/species assuming two samples from each pop by counting the number of four site types in the following order:
Ahet, Bhet, fixediff and sharedhet. Ideally, one should generalize this to arbitrary samples and # of pops.
"""

def count2X2(d):
	reslist=[]
	for bl in d.values():		
		newlist=[]
		ll=[]
		emptd = {'fixdiff':0, 'sharedhet':0, 'Ahet':0, 'Bhet':0}
		for k in bl:
			if k[0] is k[1] and k[2] is k[3]:
				t = 'fixdiff'
			elif k[0] is not k[1] and k[2] is not k[3]:
				t = 'sharedhet'
			elif k[0] is not k[1] and k[2] is k[3]:	
				t =  'Ahet'
			elif k[0] is k[1] and k[2] is not k[3]:	
				t =  'Bhet'
			newlist.append(t)
		for i in newlist:
			emptd[i] += 1
		sortk = sorted(emptd.keys())
		for key in sortk:
			ll.append(emptd[key])
		reslist.append(ll)
	return reslist


""""
Applies blockcut to a dict of chrom:coordinates and makes a new dict of blocks. The key of each block is a tuple (chrom, start, end)
"""
def main():
	block_len = int(sys.argv[4])
	ttdict = filter_dict_maxcov(bed_to_dict(sys.argv[1]),2*block_len)
	vcf_reader = vcf.Reader(filename=sys.argv[2])
	var_contig = get_chrom_set(sys.argv[3])
	samp = sys.argv[5].split(".")
	print 'no of contigs > twice the blocklength, i.e. ' + str(block_len*2)+' bases: '+str(len(ttdict.keys()))
	nndict = blockcutdict(ttdict, block_len)
	print 'no of unfiltered blocks > '+str(block_len)+' bases: '+str(len(nndict.keys()))
	ndictpostfilt = filter_dict_maxspan(nndict,2*block_len)
	print 'no of filtered blocks > '+str(block_len)+' bases + no longer than '+ str(block_len*2)+' bases: '+str(len(ndictpostfilt.keys()))
	testkeys = list(ndictpostfilt.keys())[0:10000]
	testdict = {}
	for i in testkeys:
		testdict.update({i:ndictpostfilt.get(i)})
	fetched_dict = get_vcf_contigs(vcf_reader, testdict, var_contig, min_GQ=0, min_DP=0, samples = samp)
#	fetched_dict = get_vcf_contigs(vcf_reader, ndictpostfilt, var_contig, min_GQ=0, min_DP=0, samples = samp)
#	print('no of filtered fetched blocks > 0.5kb: '+ str(len(fetched_dict.keys())))
	dict_DeNoned=filter_dict_ByNoneCount(fetched_dict,5) #could make this a sys.argv
	print('no of filtered fetched blocks with < 5 Nones: '+str(len(dict_DeNoned.keys())))
	dict_Cleaned=change_dict_RemoveMultiSNoneS(dict_DeNoned)
	print('no of filtered fetched blocks > 0.5kb with < 5 Nones + multiallelic or None sites converted to invariant: '+str(len(dict_Cleaned.keys())))
#	This is a very ugly line of code that counts up sites with 4 or more site types, which would correspond to a 4GT. At the moment we don't do anything with this info, just spit out how many non-violated blocks we have. This could also be a function if required and we could do an extra filtering step to remove these blocks.
	count4gt=[sum(1 for x in [len(set(zip(site[0],site[1]))) for site in itertools.combinations(block, 2)] if x>=4) for block in dict_Cleaned.values()]
	print('no of blocks with no 4 gamete violations '+str(count4gt.count(0)))
	np.save(sys.argv[5]+'.'+sys.argv[4]+'.npy',dict_Cleaned)
#	slist = dict_Cleaned.values()
#	out=str(slist).replace("[","{").replace("]","}").replace(")","}").replace("(","{")
#	print out
#	print 'total number of muts: '+str(count)
#	twobytwo=count2X2(dict_Cleaned)
#	out=str(twobytwo).replace("[","{").replace("]","}")
#	print out
if __name__== "__main__":
	main()