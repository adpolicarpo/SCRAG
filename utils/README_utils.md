# utils folder

**********************************
# compare_sero.py

Script to compare 2 ou 3 sets of results, obtaining the number of genes in common and the number of genes not shared (exclusive from that group)

Usage:

Open the script, and in the main() function, define path1 and path2 (path to result files), according to the examples given
Use path 3 only if you want to compare 3 datasets. In this cases you have to uncomment (remove # from) lines 34, 39, and 53 to 55.
Run script in command line: go to the folder with the script, and type "python compare_sero.py".

*******************************
# downloadFTP.py

Script to download genomes from GenBank FTP

Usage:

Open the script, and replace 'Streptococcus_pneumoniae' with the name of species wanted, written as in GenBank FTP
Run script in command line: go to the folder with the script, and type "python downloadFTP.py". Choose the extension wanted.

*******************************
# runProdigal.py

Script to run prodigal to obtain CDSs (.faa and .ffn files)

Usage:

Run script in command line: go to the folder with the script, and type "runProdigal.py 'path to contigs files'".

*******************************
# verify_clusters.py

Script to verify if a sequence is only in a cluster/CVAP

Usage:

Open the script, and in the main() function, define path_core (path to results files of core genome) and path_accessory (path to results files of accessory genome), according to the examples given
If you are using just core or just accessory genome, comment (put a # before) the rows about core/accessory genome (that one that's not being used)
Run script in command line: go to the folder with the script, and type "python verify_clusters.py".

