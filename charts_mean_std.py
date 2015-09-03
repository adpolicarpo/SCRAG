##########################################################################
##########################################################################
#################               SCRAG                    #################
#################    Strict CoRe and Accessory Genome    #################
#################        Author: Adriana Policarpo       #################
##########################################################################
##########################################################################


## Plot a chart with minimum (x), means (y) and std (color scale) of each gene, calculated using the distances matrix 


import os
import sys
import time
import glob
import numpy
import matplotlib.pyplot as plt


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

##################################################

def plot(minims, means, stds):

	# to find value next to minor of minimums, and set it as the minimum of axes
	n = []
	n2 = []
	i = 0
	while i < 1:
		n.append(round(i, 2))
		i += 0.1

	for k in n:
		if k < min(minims):
			n2.append(k)
		
	ax_min = max(n2)
		
	#plt.plot(minims, means, 'co')

	plt.scatter(x = minims, y = means, c = stds, cmap = 'jet', marker = 'o')
	plt.axis([ax_min, 1, ax_min, 1])
	cbar = plt.colorbar()
	cbar.set_label('standard deviation', rotation = 270)
	#plt.grid(axis='both')
	plt.xlabel('Min. % similarity')
	plt.ylabel('Mean')
	#plt.title('Relation between mean, standard deviation and minimum % similarity', color='r', size=12)

	figname = "min_mean_std.png"

	#plt.show()

	return figname

###################################################


def main():

	print "Starting Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")

	if analysis_type == '1':
		directory = "../Results/results_core_" + ident + os.sep + ident + "_" + dif + os.sep #directory where are  genes (results) files  ||  ident and dif = args
	elif analysis_type == '2':
		directory = "../Results/results_accessory_" + ident + os.sep + ident + "_" + dif + os.sep #directory where are the accessory genes (results) files  ||  ident and dif = args
	else:
		print "Second argument is invalid: choose '1' for core genome analysis or '2' for accessory genome analysis"


	dist_files = glob.glob(os.path.join(directory,"*.dst")) #we want only the files with the distance matrix

	gene_stats = {}

	minims = []
	means = []
	stds = []
	for dist in dist_files:
		dm = open (dist, 'r')
		lines = dm.readlines()  #retrieves a list with the lines of the file
		dm.close()
	
		nr_seq = lines[0].strip()  #number of sequences in the file
		print "File:", dist.split(os.sep)[-1], "--> corresponding to", nr_seq, "sequences"
	
		matrix = []
		i = 1
	
		line2 = []
		while i < len(lines):
			line = lines[i].strip('\n')
	
			if "gene" in line:
				if line2 != []:
					#print line2		
					matrix.append(line2)
				line2 = []
				obj = line.split(' ')
				for n in obj:
					if "." in n:
						line2.append(n)
			else:
				obj = line.split(' ')
				for n in obj:
					if "." in n:
						line2.append(n)
		
			i += 1
			if i == len(lines):
				matrix.append(line2)
			
		values = []
		j = 1  #in the first line of matrix, the value on the diagonal is 0
		while j < len(matrix):  #len(matrix) = nr_seq (number of sequences)
			line = matrix[j]
			k = 0 
			while k < len(line) and k < j:
				value = line[k]
				values.append(1 - float(value))
				k += 1
			j += 1
		
		
		minim = min(values)
		mean = numpy.mean(values)
		std = numpy.std(values)
	
		key = dist
		value = minim #, mean, std
		gene_stats[key] = value
	
		minims.append(minim)
		means.append(mean)
		stds.append(std)
	
	
		print "minimum value:", minim
		print "mean:", mean
		print "standard deviation:", std, "\n"



	###########################
	#ploting
	figname = plot(minims, means, stds)
	plt.savefig(directory + figname)
	print "Figure saved at", directory, "as", figname	
	##########################

	c = 0
	val = float(ident)/100
	print "\n\nGenes with less than " + ident + "% min similarity:\n"
	for minim in minims:
		if minim < val:
			for k, v in gene_stats.items():
				if v == minim:
					print k.split("_uniques")[0].split(os.sep)[-1]
					c += 1
					
	print "\nNumber of genes with less than " + ident + "% min similarity:", c

	len_genes = len(glob.glob(os.path.join(directory,"*uniques.fasta")))
	print "\nFiles (genes) in results data set:", len_genes
	print "After exclusion by looking for the threshold in the chart:", len_genes - c


	print "Finished Script at:", time.strftime("%H:%M:%S-%d/%m/%Y")


#######################################################################
#######################################################################


if __name__ == "__main__":
	main()


