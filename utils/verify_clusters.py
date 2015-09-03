## script to verify if a sequence is only in a cluster/CVAP

import os
import sys
import time
import glob


################################################################
##function that returns a list of gene tags in a file

def readgen (filename): 
	fn = open(filename, 'r')
	gnames = []
	i = 1
	lines = fn.readlines()
	fn.close()
	while i < len(lines):
		gnames.append(lines[i].strip())
		i += 2
	return sorted(gnames)
	
################################################################	



def main():

	print "Starting Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")

	# verify clusters in core genome, accessory genome or both
	path_core = ""	#path to results files of core -- e.g: "../Results/results_core_80/80_20/"
	path_accessory = ""	#path to results files of accessory -- e.g: "../Results/results_accessory_80/80_20/"

	#glob returns a list (even with only one element)
	results_files_core = glob.glob(os.path.join(path_core,"*.fasta"))
	results_files_accessory = glob.glob(os.path.join(path_accessory,"*.fasta"))
	
	seq_files = []
	for seq_file in results_files_core:
		if "_uniques" not in seq_file:
			seq_files.append(seq_file)
	for seq_file in results_files_accessory:
		if "_uniques" not in seq_file:
			seq_files.append(seq_file)

	#list of clusters
	genes_list = []
	for f in seq_files:
		gene_names = readgen(f)
		genes_list.append(gene_names)
	
	print "number of clusters:", len(genes_list)

	count = 0
	genes_list2 = []
	for cluster in genes_list:
		for seq in cluster:
			if seq not in genes_list2:
				genes_list2.append(seq)
			else:
				print "repeated gene: " + seq
				count += 1

	print  "Number of repeated sequences:", count
	#print len(genes_list2)


	print "Finished Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")


if __name__ == "__main__":
	main()
