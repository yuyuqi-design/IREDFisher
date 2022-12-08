#Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com
#################################################################################################################################
#function: download sequence from NCBI database or Uniprot database.
#usage: python3 fetch_sequence.py code
#example: python3 get_fasta.py P95480
#################################################################################################################################

import sys
from Bio import Entrez
from urllib import request


def fetch_sequence(code):

    try:
    #downloading the fasta individually from NCBI website using Entrez module from biopython library
        filename = str(code).strip() + '.fasta'
        with open(filename, 'w') as fasta:
            Entrez.email = "Your.Name.Here@example.org"
            handle = Entrez.efetch(db="sequences", id=str(code).strip(), rettype="fasta", retmode="text")
            fasta.write(handle.read())
    #if there is no such code in NCBI database, then download from Uniprot database
    except:
        request.urlretrieve('https://www.uniprot.org/uniprot/' + str(code).strip() + '.fasta', str(code).strip() + '.fasta')


