##########################################################################
##########################################################################
#################               SCRAG                    #################
#################    Strict CoRe and Accessory Genome    #################
#################        Author: Adriana Policarpo       #################
##########################################################################
##########################################################################


## Compares all the genomes (.faa) concatenated in a file (.fasta) with the DB with ALL genomes, by doing BLAST
## Accept all aligns with N hits (core genome analisys) or less than N (accessory genome analisys). N=number of genomes 


import sys
import os
import glob
import time
from collections import Counter
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Seq import Seq


if len(sys.argv) != 5:
	sys.stderr.write("Invalid number of arguments: 5 needed, " + str(len(sys.argv)) + " given\n")
else:
	# sys.argv[0] --> name of the script to run
	analysis_type = sys.argv[1]
	ident = sys.argv[2]
	dif = sys.argv[3]
	infile_directory = sys.argv[4]  #folder with the file to analyze (concatenated file with genomes)


IDENTITY_THRESHOLD = float(ident)/100


## Function to create BLAST databases
def Create_blastDB (input_file):
	p = input_file.split(os.sep)[-1].split('.') #the first split is to exclude de names of the folders in the directory of the file
	i = 0
	name = ''
	while p[i] != p[-1]:    #to exclude .fasta (included in the file name) of db name 
		name = name + p[i]
		i += 1
	filenames = name + "_db"
	directory = '../Results/blastDBs/' #define directory to save the database files
	if not os.path.exists(directory):
		os.mkdir(directory) #make new directory, if it not exists
	dbname = directory + filenames
	if not os.path.isfile(dbname + ".pin"):
		try:
			os.system("makeblastdb -in " + input_file + " -out " + dbname + " -dbtype prot -logfile " + dbname + "_blast.log" )
			print "Database created."
		except:
			print "Database not created."
	else:
		print "Database already existent. Using that database..."
	return dbname



## Function that makes a query to a blast database
def Query_blastDB (filename, database, results):
	if not os.path.isfile(results):
		try:
			print "Running BLAST..."
			blastp_cline = NcbiblastpCommandline(query = filename, db = database, evalue = 0.001, outfmt = 5, out = results)
			stdout, stderr = blastp_cline()
		except:
			print "Error Running BLAST"
	else:
		print "BLAST results file aready existent. Using that file..."
	result_handle = open(results)
	blast_records = NCBIXML.parse(result_handle)
	return blast_records


def Read_blast_alignments(rec, aligns, genes):
	count_idents = 0
	fnames = []
	group_list = []
	query_scores = []
	ref_score = 0
	for alignment in rec.alignments:
		hsp = alignment.hsps[0]
		max_score = hsp.score

		#hsp is that with the best score
		for match in alignment.hsps:
			if match.score > max_score:
				max_score = match.score
				hsp = match
						
		#identity = hsp.score / alignment.length
		identity = float(hsp.identities) / alignment.length
			
		if identity >= IDENTITY_THRESHOLD:
			count_idents += 1
					
		#### rstrip("[") --> Only works if the name of the strain is between "[]" ####
		origin = alignment.hit_def.split("[")[1].rstrip("]") #origin (strain) of the sequence. 
		fnames.append(origin)

		#list of genes definitions
		group_list.append(alignment.hit_def.strip())


		#reference score to get Blast Score Ratio.
		#The ref score is the score of the alignment of a sequence against itself
		#in each rec there are 1 ref score and "NFILES" query scores (one of them is the same that ref score)
		if alignment.hit_def.strip() == rec.query.strip():
			ref_score = hsp.score
		query_score = hsp.score
		query_scores.append(query_score)


	#calculate the Blast Score Ratio (BSR)
	#we get "NFILES" BSRs for each blast record, one of them correspondent to the query agains itself (=1)
	BSR_cutoff = True
	for query_score in query_scores:
		BSR = query_score/ref_score
		if BSR < 0.6:
			BSR_cutoff = False
			break


	# We will now determine if all aligns in the record are from the different "sources" (strains)
	# If "Counter" doesn't work with the python's version installed, use next block instead of this one	
	##########################
	stop = False
	c = Counter(fnames)	#counts all the ocurrences of each element in the list 
	for value in c.values():
		if value != 1:	#we want only one of each (means that all the strains are different)
			stop = True
			break
	##########################
	# If "Counter" doesn't work with the python's version installed, use this block instead of the previous one
	##########################
	# stop = False
	# 	i = 0
	# 	while i < len(fnames)-1 and stop == False:
	# 		j = i + 1
	# 		while j < len(fnames) and stop == False:
	# 			if fnames[i] == fnames[j]:
	# 				stop = True
	# 			j += 1
	# 		i += 1
	##########################

	
	if count_idents == len(rec.alignments) and stop == False and sorted(group_list) not in genes and BSR_cutoff == True:
		aligns.append(group_list)
		genes.append(sorted(group_list))

		# each record is considered only if:
		# - have "NFILES" alignments (core) or less than NFILES alignments (accessory)
		# - all the alignments on it have the given threshold
		# - all the sequences(alignments) are from diferent sources (strains)
		# - each record is unique (not repeated), and
		# - BSRs values in the set of alignments is >= 0.6

	return aligns, genes


##Function to read blast records
def Read_blast_rec (blast_records, directory):
	NFILES = len(os.listdir(directory)) - 1  #directory is the folder with .faa files and the concatenated .fasta
	if analysis_type == '1':
		aligns = []
		genes = []
		for rec in blast_records:
			#print len(rec.alignments)
			if len(rec.alignments) == NFILES: #only the records with lenght = NFILES ("NFILES" alignments) are considered
				aligns, genes = Read_blast_alignments(rec, aligns, genes)
		print "Found " + str(len(aligns)) + " genes with " + str(NFILES) + " different matches, with > " + str(ident) + "%" + " similarity and BSR >/= 0.6\n"


	elif analysis_type == '2':
		aligns = []
		genes = []
		for rec in blast_records:
			#print len(rec.alignments)
			if len(rec.alignments) < NFILES: #only the records with lenght < NFILES ("NFILES" alignments) are considered
				aligns, genes = Read_blast_alignments(rec, aligns, genes)
		print "Found " + str(len(aligns)) + " genes with less than " + str(NFILES) + " different matches, with > " + str(ident) + "%" + " similarity and BSR >/= 0.6\n"
	
	return aligns



## function to select only sequences with less than x% size difference (in each group)
def exclude_by_size(groups_of_genes, dict_tags_seqs, size_dif):
	count = 0
	groups_filtered = []
	for group in groups_of_genes:
		sizes = []
		for gen in group:
			size = len(dict_tags_seqs[gen])
			sizes.append(size)
		minim = min(sizes)
		maxim = max(sizes)
		ratio = float(minim)/maxim

		if ratio >= size_dif:
			groups_filtered.append(group)
		else:
			count += 1

	return groups_filtered, count


###################################################################
###################################################################



def main():

	print "Starting Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")

	print "\nIdentity threshold:", IDENTITY_THRESHOLD, "(" + ident + "%)"
	print "Maximum size difference between sequences allowed: " + dif + "%\n"


	#infile_directory is an argument

	#only the concatenated file is .fasta, the remaining are .faa. glob returns a list (even with only one element)
	infile = glob.glob(os.path.join(infile_directory,"*.fasta"))[0]
	
	# directory to save the results
	directory = "../Results/" #all the results will be saved in this folder
	if not os.path.exists(directory):
		os.mkdir(directory)


	if analysis_type == '1':
		newdir = directory + "results_core_" + ident + "/"  #folder to save the results of core genome analysis
		if not os.path.exists(newdir):
			os.mkdir(newdir)
		print "Performing a core genes analysis\n"
	elif analysis_type == '2':
		newdir = directory + "results_accessory_" + ident + "/"  #folder to save the results of accessory genome analysis
		if not os.path.exists(newdir):
			os.mkdir(newdir)
		print "Performing an accessory genes analysis\n"
	else:
	 	print "Second argument is invalid: choose '1' for core genome analysis or '2' for accessory genome analysis"


	#Directory to save the results before the size exclusion
	no_excsize = newdir + "no_size_exclusion" + os.sep
	if not os.path.exists(no_excsize):
		os.mkdir(no_excsize)

	#Create genes dictionary from file with all genomes, to recover original sequences
	fl = open(infile, "rt")
	lines = fl.readlines()
	allgenes = {} # empty dictionary
	defin = ''
	seq = ''
	for line in lines:
		if ">" in line:
			if seq and defin != '':
				allgenes[defin] = seq.rstrip("*")
			defin = line.strip().lstrip(">")
			seq = ''
		else:
			seq += line.strip()
	allgenes[defin] = seq.rstrip("*")
	fl.close()


	############## only do this block if it's the first analysis for this %identity, otherwise use existent files ############
	
	if len(glob.glob(os.path.join(no_excsize,"*.fasta"))) == 0:

		#the file to query the database is the same that creates it
		#Create a blast db
		db_file = infile
		db = Create_blastDB(db_file)

		# Query a blast db and parsing the results
		query_file = infile
		blast_results_file = db + "_blast_results.xml"
		blast_recs = Query_blastDB (query_file, db, blast_results_file)

		#Filtering the blast results according to parameters given
		aligns = Read_blast_rec (blast_recs, infile_directory)	


		if aligns != '':

			#save results in file (before exclusion by size)
			if analysis_type == '1':
				rname = "coregene"
			elif analysis_type == '2':
				rname = "accessorygene"
			i = 0
			for group in aligns:
				i += 1
				result = open(no_excsize + rname + str(i) + ".fasta", 'w') #open a file to write
				result.write("Found a gene with " + str(len(group)) + " match(es)\n") #len(group) is the number of sequences in an alignment
				for defin in group:
					seq = Seq(allgenes[defin])
					result.write (">" + defin + "\n") #definition of matching sequence
					result.write (str(seq) + "\n") #matching sequence, (without gaps, recovered from de dict with all genes)
				result.close()

			print "Saved in files before exclusion by size difference"

	else:
		print "Files with more than " + ident + "% similarity already existent. Proceding now to exclusion by size."
		print len(glob.glob(os.path.join(no_excsize,"*.fasta"))), "genes found."


		###################################################################################################################


	########## Exclude genes by size difference ###########
	files_genes = glob.glob(os.path.join(no_excsize,"*.fasta"))

	read_genes = []
	for f in files_genes:
		group = []
		result = open(f, "rt")
		lines = result.readlines()
		for line in lines:
			if ">" in line:
				group.append(line.strip().lstrip(">"))
		read_genes.append(group)
		result.close()


	SIZE_DIF = 1 - (float(dif)/100)  #sequences must not have more than 1-SIZE_DIF (SIZE_DIF%) difference in length	
	groups_genes, nr_excluded = exclude_by_size(read_genes, allgenes, SIZE_DIF)

	print nr_excluded, "genes excluded by size variation > " + dif + "%"
	print "Number of genes with less than " + dif + "% variation in size:", len(groups_genes), "\n"
	################

	#save results in file (after exclude by size difference)
	dir_excsize = newdir + ident + "_" + dif + "/"
	if not os.path.exists(dir_excsize):
		os.mkdir(dir_excsize)

	if analysis_type == '1':
		rname = "coregene"
	elif analysis_type == '2':
		rname = "accessorygene"
	i = 0
	for group in groups_genes:
		i += 1
		result = open(dir_excsize + rname + str(i) + ".fasta", 'w') #open a file to write
		result.write("Found a gene with " + str(len(group)) + " match(es)\n") #len(group) is the number of sequences in an alignment
		for defin in group:
			seq = Seq(allgenes[defin])
			result.write (">" + defin + "\n") #definition of matching sequence
			result.write (str(seq) + "\n") #matching sequence, (without gaps, recovered from de dict with all genes)
		result.close()
	

	print "Finished Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")


#######################################################################
#######################################################################

if __name__ == "__main__":
	main()

