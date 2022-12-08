import subprocess as sp
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pandas as pd
import sys


def parse_blast(protein='', file='blast.xml', single_organism=False, save_csv=True):
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

    for blast_record in blast_records:
        #print(len(blast_records.alignments))
        #every blast_records has three attributes: descriptions, alignments and multiple_alignment
        for alignment in blast_record.alignments:
            #alignment attributes of blast_records: title, length and hsps (a list containing all seqeunces)
            #example for alignment.title: ref|WP_006750377.1| tryptophan 7-halogenase [Burkholderia ambifaria]
            #we can extract the Organism, Accession ID, Descripiton from title attributes.
            #for example, the orsm will be Burkholderia_ambifaria
            orsm= '_'.join(alignment.title.strip().split(']')[-2].split('[')[-1].split())
            #print(orsm)
            #acesion is WP_006750377
            acesion = alignment.title.strip().split('|')[1].split('.')[0]
            #print(acesion)
            #descrip is tryptophan 7-halogenase
            descrip = alignment.title.strip().split('|')[2].split('[')[0]
            #add each item to the list
            Organism.append(orsm)
            Accession.append(acesion)
            Description.append(descrip)

            for hsp in alignment.hsps:
                print(hsp.num_alignments)
                identity = hsp.identities/hsp.align_length
                #if (identity > 0.3) and (hsp.expect < 0.01):
                E.append(hsp.expect)
                Identity.append(identity)
    #print(len(E))
    #print(len(Organism))
    Accession = pd.Series(Accession, name='Accession')
    Description = pd.Series(Description, name='Description')
    Organism = pd.Series(Organism, name='Organism')
    E = pd.Series(E,name='E-value')
    Identity = pd.Series(Identity,name='Identity')
    output = pd.concat([Accession,Description, Organism,E, Identity], axis = 1)
    #print(output)

    # get those sequences with evalue greater than 0.01 and identity greater than 0.3
    #only keep those sequences whose E value < 10 and sequence identity over 30% 
    output_filter_e_identity=output[(output['E-value']< 10) & (output['Identity'] > 0.3)]
    #if the user specify the protein name
    if protein:
        df_namefilt = output_filter_e_identity[output_filter_e_identity['Description'].str.contains(protein)]
    else:
        df_namefilt = output_filter_e_identity
    #if user choose to select one sequence from a single organism
    df_unqOrgm = df_namefilt.drop_duplicates(subset='Organism', keep='first')

    #Header: Accession,Description, Organism,E, Identity
    if save_csv == True:
        df_unqOrgm.to_csv('blast.csv',index=False)
    else:
        return df_unqOrgm
    result_handle.close() 