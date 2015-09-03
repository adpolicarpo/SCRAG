##########################################################################
##########################################################################
#################               SCRAG                    #################
#################    Strict CoRe and Accessory Genome    #################
#################        Author: Adriana Policarpo       #################
##########################################################################
##########################################################################


##   Read gene files and obtain unique sequences with a new tag
## (recovering DNA sequences, intead of using aminoacid sequences) 


import os
import sys
import glob
import time
from multiprocessing import Pool


if len(sys.argv) != 6:
	sys.stderr.write("Invalid number of arguments: 6 needed, " + str(len(sys.argv)) + " given\n")
	analysis_type = ''
	ident = ''
	dif = ''
else:
	# sys.argv[0] --> name of the script to run
	prot_directory = sys.argv[1]
	nucl_directory = sys.argv[2]
	analysis_type = sys.argv[3]
	ident = sys.argv[4]
	dif = sys.argv[5]


def GetDNAseqs(protseqsfile, DNAseqsfile):
	convert_prot_nucl = {} # empty dictionary
	protdefins = []
	nuclseqs =[]
	prot = open(protseqsfile, "rt")
	nucl = open(DNAseqsfile, "rt")
	prots = prot.readlines()
	nucls = nucl.readlines()

	for line in prots:
		if ">" in line:
			protdefins.append(line.strip().lstrip(">"))
	prot.close()


	nuclseq = ''
	for line in nucls:
		if ">" in line:
			if nuclseq != '':
				nuclseqs.append(nuclseq)
				nuclseq = ''
		else:
			nuclseq += line.strip()
	nuclseqs.append(nuclseq)
	nucl.close()

	#print len(protdefins), len(nuclseqs)
	i = 0
	for defin in protdefins:
		convert_prot_nucl[defin] = nuclseqs[i]
		i += 1

	return convert_prot_nucl


def GetUniqueSeqs(seqs_file):
	read = open(seqs_file, 'r')
 	lines = read.readlines()
 	read.close()
 	seq_uniques = {}
 	i = 2 #line with the first sequence
 	while i < len(lines):
 		seq_uniques[lines[i].strip()] = lines[i-1].strip().strip(">") #no duplicate keys in a dictionary, so, using seqs as keys, we will get only the different sequences
 		i += 2  #lines with sequences
 	read.close()

 	name = seqs_file.split(".")[-2]
 	newfilename = ".." + name + "_uniques.fasta"  #.. is because of the dir name and  [-2] to exclude .fasta
 	newf = open(newfilename, 'w')
	x = 1
	for prot_defin in seq_uniques.values():
		tag = ">" + name.split(os.sep)[-1] + "_" + str(x) + "\n"
		seq = translator[prot_defin]
		newf.write(tag)
 		newf.write(seq + "\n")
 		x += 1
 	newf.close()


####################################################


def main():

	print "Starting Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")

	#prot_directory --> sys.argv[1]
	#only the concatenated file is .fasta, the remaining are .faa. glob returns a list (even only one element)
	prot_file = glob.glob(os.path.join(prot_directory,"*.fasta"))[0]

	#nucl_directory --> sys.argv[2]
	#only the concatenated file is .fasta, the remaining are .faa. glob returns a list (even only one element)
	nucl_file = glob.glob(os.path.join(nucl_directory,"*.fasta"))[0]

	#to get DNA sequences correspondent to protein sequences
	global translator
	translator = GetDNAseqs(prot_file, nucl_file)

	
	if analysis_type == '1':
		directory = "../Results/results_core_" + ident + os.sep + ident + "_" + dif + os.sep #directory where are the core genes (results) files  ||  ident and dif = args
	elif analysis_type == '2':
		directory = "../Results/results_accessory_" + ident + os.sep + ident + "_" + dif + os.sep #directory where are the accessory genes (results) files  ||  ident and dif = args
	else:
		print "Second argument is invalid: choose '1' for core genome analysis or '2' for accessory genome analysis"

	files = glob.glob(os.path.join(directory,"*.fasta")) 
	seq_files = []
	for seq_file in files:
		if "_uniques" not in seq_file:
			seq_files.append(seq_file)

	
	poolJobs = Pool()
	poolJobs.map(GetUniqueSeqs, seq_files)
	poolJobs.close()

	print "Finished Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")


####################################################


if __name__ == "__main__":
	main()
