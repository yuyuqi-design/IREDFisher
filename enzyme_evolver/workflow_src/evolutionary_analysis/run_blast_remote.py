#Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com
####################################################################################################################################################################
#function: run and parse blast result                                                                                                                               
#requirement: refseq_protein.00 NCBI database
#usage: python3 run_blast.py query_sequence(fasta format) protein_name
#example: python3 run_blast.py P95480.fasta halogenase
#output: blast.xml, blast.csv (a summary file)
#####################################################################################################################################################################

import subprocess as sp
import sys
import uuid


# ###def run_blast(query, output_name='blast.xml', blastp='blastp', refseqdb='swissprot'):
    ##############################################blast locally   
    #Blast in NCBI refseq_protein database or in-house database(upload by users)
    #input: query 
    #output: a xml file 
    #-outfmt equals 5 means an xml format result
    #option1: specify the protein name (eg, halogenase)
    #option2 : one sequence one organism
    #option3: upload an in-house database (later) by converting fasta format to database format using 
    #         makeblastdb -in inhouse.fasta -parse_seqids -blastdb_version 5 -title "inhouse" -dbtype prot
    #############################################blast remotely
def run_blast_remote(query, refseqdb='refseq_protein'):
    # fasta_string = open(query).read()
    # result_handle = NCBIWWW.qblast('blastp', refseqdb, fasta_string,hitlist_size=400, results_file='blast.xml')
    command = 'blastp' + ' -query ' + query + ' -db ' + refseqdb + ' -out ' + 'blast.xml' + ' -outfmt 5 -max_hsps 1 -max_target_seqs 400 -remote'
    sp.run(command, shell=True)
    ####parse the result

