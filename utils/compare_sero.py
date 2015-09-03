## script to compare 2 ou 3 sets of results, obtaining the number of genes in common
## and the number of genes not shared (exclusive from that group)


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

	path1 = ""	#path to results 1 (less genomes) -- e.g.: "./ser1_3/Results/results_accessory_80/80_20/"
	path2 = ""	#path to results 2 (more genomes) -- e.g.: "./all76gen/Results/results_accessory_80/80_20/"
	#path3 = ""  #path to results 3 (more genomes) -- e.g.:"./outros_s3/Results/results_accessory_80/80_20/" -- only use this when comparig 3 sets of results

	#glob returns a list (even with only one element)
	results_files1 = glob.glob(os.path.join(path1,"*.fasta"))
	results_files2 = glob.glob(os.path.join(path2,"*.fasta"))
	#results_files3 = glob.glob(os.path.join(path3,"*.fasta"))	##only use this when comparig 3 sets of results

	seq_files1 = []
	for seq_file in results_files1:
		if "_uniques" not in seq_file:
			seq_files1.append(seq_file)
	seq_files2 = []
	for seq_file in results_files2:
		if "_uniques" not in seq_file:
			seq_files2.append(seq_file)
	

	##only use this block when comparig 3 sets of results
	##########
	#for seq_file in results_files3:
	#	if "_uniques" not in seq_file:
	#		seq_files2.append(seq_file)
	#########


	genes_list = []
	for f in seq_files1: 	 #we want genes in set1, but not in set2 (and set 3)
		gene_names1 = readgen(f)
		genes_list.append(gene_names1)
	print "number of clusters in set1:", len(genes_list)

	genes_list2 = []
	for f in seq_files2:	
		gene_names2 = readgen(f)
		genes_list2.append(gene_names2)
	print "number of clusters in set2:", len(genes_list2)

		
	count = 0
	count2 = 0
	for cluster in genes_list:
  		count_temp = count
 		for cluster2 in genes_list2:
 			if len(cluster2) > len(cluster):
 				dif = set(cluster2) - set(cluster)  #list with different genes
 				if len(dif) == len(cluster2) - len(cluster):
    				count += 1
    		else:
    			dif = set(cluster) - set(cluster2)  #list with different genes
 				if len(dif) == len(cluster) - len(cluster2):
    				count += 1
  		if count == count_temp:
    		count2 += 1

	print "genes in common:", count
	print "genes not in common:", count2



	print "Finished Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")



if __name__ == "__main__":
	main()

