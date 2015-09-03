#########################################################
#####  Script to download genomes from GenBank FTP  #####
#####           Author: Adriana Policarpo           #####
#########################################################


from ftplib import FTP
import os

ftp = FTP('ftp.ncbi.nih.gov')   # connect to host
ftp.login()
ftp.cwd('genomes/Bacteria/') # set directory

l = ftp.nlst()  # retrieves a list of file names

print 'Matching Streptococcus_pneumoniae...'
print 'Select the extension of the files you wish to download'
print 'e.g. .gff, .tab, .val, .asn, .fna...'

resp = raw_input()

#############

# function that returns the name of the files with the choosen extention,
# in the given directory
def files (name, ext):
    ftp.cwd(name) # set directory
    fil = ftp.nlst()   # retrieves a list of file names
    for f in fil:
        if ext in f:
            return f
        
############

# replace 'Streptococcus_pneumoniae' with the name of species wanted, written as in GenBank FTP
for item in l:
    if 'Streptococcus_pneumoniae' in item:
        try:
            filename = files (item, resp)
            filename2 = item + '_' + filename
            print filename2
            directory = './'+item #define directory to save the files
            if not os.path.exists(directory):
                os.mkdir(directory) #make new directory, if it not exists
            os.chdir(directory) #change directory
            file = open(filename2, 'wb') #open a file
            ftp.retrbinary('RETR ' + filename, file.write) #download file
            file.close()    #close the file
            ftp.cwd('..')    # returns to the last directory
            os.chdir('..')
        except:
            print 'no files with this extention'


ftp.close() # Close the connection
    
    
