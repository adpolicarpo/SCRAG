# script to run prodigal to obtain CDSs (.faa and -ffn files)

import sys
import os
import fileinput

def main():
	
	try:
		input_path = sys.argv[1]
	except:
		print "usage: path to contigs files"
	
	
	nucl_path = "./prodigal_nucl_" + input_path.split(os.sep)[1] + os.sep
	if not os.path.exists(nucl_path):
		os.mkdir(nucl_path)
	prot_path = "./prodigal_prot_" + input_path.split(os.sep)[1] + os.sep
	if not os.path.exists(prot_path):
		os.mkdir(prot_path)


	input_files = os.listdir(input_path)
	for input_file in input_files:
		contigsFasta = input_path + input_file
		nucl_file = nucl_path + input_file.split(".")[0] + ".ffn"
		prot_file = prot_path + input_file.split(".")[0] + ".faa"

		#RUN PRODIGAL
		prod_command = "prodigal " + '-i ' + contigsFasta + ' -c ' + '-m ' + '-g 11 ' + '-p single ' + '-q ' + '-d ' + nucl_file + ' -a ' + prot_file
		print prod_command
		run_prodigal = os.system(prod_command)

		strain_name  = input_file.split(".")[0]
		with open(nucl_file, "r+") as nucl:
			for line in fileinput.input(nucl_file, inplace = True):
				if ">" in line:
					newline = line.split('#')[0] + "[" + strain_name + "]"
					print newline
				else:
					print line.strip()

		with open(prot_file, "r+") as prot:
			for line in fileinput.input(prot_file, inplace = True):
				if ">" in line:
					newline = line.split('#')[0] + "[" + strain_name + "]"
					print newline
				else:
					print line.strip()


if __name__ == "__main__":
    main()
