##########################################################################
##########################################################################
#################               SCRAG                    #################
#################    Strict CoRe and Accessory Genome    #################
#################        Author: Adriana Policarpo       #################
##########################################################################
##########################################################################

## script to run all the scripts to performe core genome analysis or accessory genome analysis

import os
import sys
import time


if len(sys.argv) != 4:
	sys.stderr.write("Invalid number of arguments: 4 needed, " + str(len(sys.argv)) + " given\n")
else:
	#analysis type:
		#core genome analysis = 1, accessory genome analysis = 2
	# sys.argv[0] --> name of the script to run
	analysis_type = sys.argv[1]
	identity_threshold = sys.argv[2]
	size_dif_max = sys.argv[3]


# function to concatenate all files to analize in one file
def ContactenateFiles(files_path, output_file):
	files = sorted(os.listdir(files_path))
	with open(output_file, 'w') as outfile:
		for fname in files:
			with open(files_path + fname) as infile:
				for line in infile:
					outfile.write(line)



def main():

	print "Started to run:", time.strftime("%H:%M:%S-%d/%m/%Y")

	#####################################################################
	#Create the concatenated files (protein and DNA sequences)
	prot_directory = "../faa_files/"
	output_prot_file = prot_directory + "all_genomes_prot.fasta"
	if not os.path.isfile(output_prot_file):
		ContactenateFiles(prot_directory,output_prot_file)
		print "Concatenated file with amino acid sequences created."
	else:
		print "Concatenated file with amino acid sequences already existent."

	nucl_directory = "../ffn_files/"
	output_nucl_file = nucl_directory + "all_genomes_nucl.fasta"
	if not os.path.isfile(output_nucl_file):
		ContactenateFiles(nucl_directory, output_nucl_file)
		print "Concatenated file with nucleotide sequences created."
	else:
		print "Concatenated file with nucleotide sequences already existent."

	#####################################################################
		
	#run all the scripts of SCRAG
	os.system("python gen_analysis.py " + analysis_type + " " + identity_threshold + " " + size_dif_max + " " + prot_directory)
	os.system("python get_unique_seqs.py " + prot_directory + " "  + nucl_directory + " "  + analysis_type + " " + identity_threshold + " " + size_dif_max)
	os.system("python msa.py " + analysis_type + " " + identity_threshold + " " + size_dif_max)
	os.system("python charts_mean_std.py " + analysis_type + " " + identity_threshold + " " + size_dif_max)



	print "Ended to run:", time.strftime("%H:%M:%S-%d/%m/%Y")



if __name__ == "__main__":
	main()


