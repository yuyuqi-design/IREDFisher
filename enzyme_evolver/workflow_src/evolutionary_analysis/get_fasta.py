#Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com
###############################################################################
#This function is to download the sequence in Fasta format from NCBI database
# Input: a txt extension file containing the a list of accession numbers OR a single accession number              
# Output: all fasta files for all accession numbers
# usage: python3 get_fasta.py file.txt  OR python3 get_fasta.py P95480
###############################################################################

import sys
import subprocess as sp
import os
from Bio import Entrez
from urllib import request

def fetch_sequence(code,ofolder):
    filename = str(code).strip() + '.fasta'
    if os.path.isfile(ofolder+filename):
            pass
    else:
        #function: download sequence from NCBI database or Uniprot database.
        try:
        #downloading the fasta individually from NCBI website using Entrez module from biopython library 
            
            with open(ofolder+filename, 'w') as fasta:
                Entrez.email = "Your.Name.Here@example.org"
                handle = Entrez.efetch(db="sequences", id=str(code).strip(), rettype="fasta", retmode="text")
                fasta.write(handle.read())
        #if there is no such code in NCBI database, then download from Uniprot database
        except:
            request.urlretrieve('https://www.uniprot.org/uniprot/' + str(code).strip() + '.fasta', ofolder+str(code).strip() + '.fasta')

###############################################################################################################################

def get_fasta(file,ofolder=''):

    # check if the input is a file format.
    if str(file[-1:-5:-1][::-1]) == '.txt':
    #If it is a file with extension .txt
        with open(file) as unicodes:
            all_codes = unicodes.readlines()
            for code in all_codes:
                code = code.strip()
                fetch_sequence(code,ofolder)
    else:
        code = str(file).strip()
        fetch_sequence(code.strip(),ofolder)
def clean_fasta(file,ofolder=None):
    # check if the input is a file format.
    if str(file[-1:-5:-1][::-1]) == '.txt':
    #If it is a file with extension .txt
        with open(file) as unicodes:
            all_codes = unicodes.readlines()
            for code in all_codes:
                sp.run('rm ' + ofolder+ code.strip() +'.fasta ', shell=True)
             #   fetch_sequence(code)
    else:
        code = str(file).strip()
        sp.run('rm ' + ofolder +code +'.fasta',shell=True)

    
        
################################################################################################################################
if __name__ == '__main__':
    file = sys.argv[1]
    get_fasta(file)
