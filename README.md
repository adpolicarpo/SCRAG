# SCRAG
Strict CoRe and Accessory Genome -  scripts developed during the master's thesis


# Introduction

SCRAG (Strict CoRe and Accessory Genome) is a bioinformatic tool, allowing the comparison of several genomes simultaneously, obtaining the core and accessory genome in a conservative way. SCRAG is based on the sequence comparison process using the BLAST algorithm, whose results are then filtered by various parameters of which the user can define the percentage of identity and the percentage of maximum size difference allowed between sequences of a set of alleles encoding the same locus (CVAPs - "Conjuntos de Variantes Alélicas Possíveis").

SCRAG also performs a multiple sequence alignment using MUSCLE and obtains a distance matrix using ClustalW for each CVAP. Using the inverse of the distances (1 - distance), it calculates the similarity valus, and the respective minimum, mean and standard deviation between sequences of the same CVAP. Using this statistics, it builds a dot plot of the mean in function of the minimum, and standard deviation as color scale. Each point represents a CVAP.


# Structure of the files to use with SCRAG

    >Tag of the sequence
    Sequence
    
e.g.:

    >gi|307126169|ref|YP_003878200.1| hypothetical protein SP670_0018 [Streptococcus pneumoniae 670-6B]
    MEQYTKRASALIITFNRGVISQDEFIEEFTKLREKVVKAVSEAKNDKL

* Files with the extension .faa (proteomes) and .ffn(genomes) should be used.
* Lines with tags should have ">" before the description of the sequence
* .faa files and .ffn files should correspond to the same proteomes/genomes and have the sequences corresponding in the same order.
* .faa files should have the origin (strain) of the sequence between brackets "[ ]" (e.g.: [Streptococcus pneumoniae 670-6B])


# Usage

1. Download SCRAG
2. In the folder "faa_files" put the .faa files to analyze, and in the folder "ffn_files" put the .ffn files. 
3. In command line, go to the directory with the scripts of the program, named "scripts" (where the "run_all.py" file is).
4. Type 'python run_all.py "type_of_analysis" "percentage of identity" "percentage of size difference"' (e.g.: 'python run_all.py 1 80 20')

    * "type_of_analysis" - 1 = core genome, 2 = accessory genome
    * "percentage of identity" - percentage of identity between sequences of a CVAP
    * "percentage of size difference" - percentage of maximum size difference allowed between sequences of a CVAP

