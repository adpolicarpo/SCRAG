##########################################################################
##########################################################################
#################               SCRAG                    #################
#################    Strict CoRe and Accessory Genome    #################
#################        Author: Adriana Policarpo       #################
##########################################################################
##########################################################################


## Performs MSA using MUSCLE and then calculates the distances matrix using ClustalW


import os
import sys
import time
import glob
from multiprocessing import Pool
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalwCommandline


if len(sys.argv) != 4:
	sys.stderr.write("Invalid number of arguments: 4 needed, " + str(len(sys.argv)) + " given\n")
	analysis_type = ''
	ident = ''
	dif = ''
else:
	# sys.argv[0] --> name of the script to run
	analysis_type = sys.argv[1]
	ident = sys.argv[2]
	dif = sys.argv[3]


###################################################

#To perform MSA using muscle
def MuscleMSA(gene_file):
	inp = gene_file
	outp = ".." + gene_file.split('.')[-2] + "_align.aln"  #.. is because of the dir name and  [-2] to exclude .fasta
	muscle_cline = MuscleCommandline(input = inp, out = outp, clwstrict = True)
	#print muscle_cline
	stdout, stderr = muscle_cline()


#To perform a tree, obtaining also distance matrices
def Clustal_dist_matrix(align_file):
	clw_cline = ClustalwCommandline("clustalw", infile = align_file, tree = True, outputtree = "dist") #outputtree = "dist" --> gives us the distances matrix
	#print clw_cline
	stdout, stderr = clw_cline()


def main():

	print "Starting Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")

	if analysis_type == '1':
		directory = "../Results/results_core_" + ident + os.sep + ident + "_" + dif + os.sep #directory where are  genes (results) files  ||  ident and dif = args
	elif analysis_type == '2':
		directory = "../Results/results_accessory_" + ident + os.sep + ident + "_" + dif + os.sep #directory where are the accessory genes (results) files  ||  ident and dif = args
	else:
		print "Second argument is invalid: choose '1' for core genome analysis or '2' for accessory genome analysis"

	print "Performing MSA with MUSCLE..."
	gene_files = glob.glob(os.path.join(directory,"*uniques.fasta")) #we want only the files with unique sequences (not with sequences repeated)
	poolJobs = Pool()
	poolJobs.map(MuscleMSA, gene_files)
	poolJobs.close()
	poolJobs.join()
	print "Obtaining distances matrices using ClustalW..."
	align_files = glob.glob(os.path.join(directory,"*.aln")) # --> now we wnat only the new files generated with muscle
	poolJobs = Pool()
	poolJobs.map(Clustal_dist_matrix, align_files)
	poolJobs.close()

	print "Finished Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")


#######################################################################
#######################################################################


if __name__ == "__main__":
	main()
