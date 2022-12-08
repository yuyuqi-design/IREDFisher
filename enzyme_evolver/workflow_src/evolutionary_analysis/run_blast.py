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
from Bio.Blast import NCBIXML
import pandas as pd
import sys


def run_blast(query, output_name='blast.xml', blastp='blastp', refseqdb='swissprot'):
    #Blast in NCBI reference sequence database
    #put refseq_protein.00 in current directory
    #input: query accession number, output: a xml file containing blast result which can be opened by a browser
    #-outfmt equals 5 means an xml format result
    command = blastp + ' -query ' + query + ' -db ' + refseqdb + ' -out ' + output_name + ' -outfmt 5 -max_hsps 1 -remote'
    sp.run(command, shell=True)
#######################################################################################################################################################################
def parse_blast(protein, file='blast.xml', save_csv=True):
    #Parse the blast result
    #input: blast xml result; the name of the enzyme. for example: halogenase
    # output: a csv file consisting of the Accession ID number, Description (protein name), Organism, E-value, sequence identity between query and the hit sequence
    result_handle = open(file)
    blast_records = NCBIXML.parse(result_handle)
    #blast_records contains all high scoring sequences
    blast_records = list(blast_records)
    # print(len(blast_records))
    #create lists for Accession ID number, Description (protein name), Organism, E-value, sequence identity between query and the hit sequence.
    #all sequences have e-value<0.01, sequence identity > 0.3 and one organism only keep one protein.
    Accession = []
    Description = []
    Organism = []
    E = []
    Identity = []

    for blast_records in blast_records:
    #     print(len(blast_records.alignments))
        #every blast_records has three attributes: descriptions, alignments and multiple_alignment
        for alignment in blast_records.alignments:
            #alignment attributes of blast_records: title, length and hsps (a list containing all seqeunces)
            #example for alignment.title: ref|WP_006750377.1| tryptophan 7-halogenase [Burkholderia ambifaria]
            #we can extract the Organism, Accession ID, Descripiton from title attributes.
            #for example, the orsm will be Burkholderia_ambifaria
            orsm= '_'.join(alignment.title.strip().split(']')[-2].split('[')[-1].split())
            #acesion is WP_006750377
            acesion = alignment.title.strip().split('|')[1].split('.')[0]
            #descrip is tryptophan 7-halogenase
            descrip = alignment.title.strip().split('|')[2].split('[')[0]
            #add each item to the list
            Organism.append(orsm)
            Accession.append(acesion)
            Description.append(descrip)

            for hsp in alignment.hsps:
    #             print(hsp.num_alignments)
                identity = hsp.identities/hsp.align_length
                #if (identity > 0.3) and (hsp.expect < 0.01):
                E.append(hsp.expect)
                Identity.append(identity)
    # print(len(E))
    # print(len(Organism))
    Accession = pd.Series(Accession, name='Accession')
    Description = pd.Series(Description, name='Description')
    Organism = pd.Series(Organism, name='Organism')
    E = pd.Series(E,name='E-value')
    Identity = pd.Series(Identity,name='Identity')
    output = pd.concat([Accession,Description, Organism,E, Identity], axis = 1)

    # get those sequences with evalue greater than 0.01 and identity greater than 0.3
    #only keep those sequences whose E value < 10 and sequence identity over 30% 
    output_filter_e_identity=output[(output['E-value']< 10) & (output['Identity'] > 0.3)]
    #extract the description is halogenase
    df_halogenase = output_filter_e_identity[output_filter_e_identity['Description'].str.contains(protein)]
    #select one which has biggest sequence identity for each organism
    df_unqOrgm = df_halogenase.drop_duplicates(subset='Organism', keep='first')

    #Header: Accession,Description, Organism,E, Identity
    if save_csv == True:
        df_unqOrgm.to_csv('blast.csv',index=False)
    else:
        return df_unqOrgm
    result_handle.close() 
################################################################################################################################

query = sys.argv[1]
# protein_name = sys.argv[2]
run_blast(query)
# parse_blast(protein_name)
