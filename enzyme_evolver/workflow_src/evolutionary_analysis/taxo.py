#Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com
####################################################################################################################################################################
#function: Taxonomy annotation based on NCBI taxonomy database      
#Input: hit_summary.csv from hit_summary.py                                                                                                                 
#Usage: python3 tax.py hit_summary.csv
#output: hit_summary_taxo.csv (header: Tax ID,Organism,Tax_kingdom,Tax_group,Accession,Description,E-value,Identity,Sequence,Tax_blastName)
#####################################################################################################################################################################
import pandas as pd
import sys
import csv
import subprocess as sp

def taxo(df=None,file='', save_csv=True):
    if file == '':
        df_sum = df
    else:
        #get the organism name from selected seqeunces
        df_sum = pd.read_csv(file)

    organism = df_sum['Organism'].str.strip().str.replace('_', ' ')
    #print(organism)
    #NCBI taxonomy database: taxID, organism,taxonomy, nan
    df_taxoNCBI = pd.read_csv('/home/g02808yy/data/database/fullnamelineage.dmp', delimiter='|', names = ['TaxID','Organism','Tax','nan'])
    #remove tab from organism and tax collum
    df_taxoNCBI['Organism'] = df_taxoNCBI['Organism'].str.strip()
    df_taxoNCBI['Tax_kingdom'] = df_taxoNCBI['Tax'].str.strip().str.split(';',expand = True)[1]
    df_taxoNCBI['Tax_group'] = df_taxoNCBI['Tax'].str.strip().str.split(';',expand = True)[4]
    #extract sequences organisms from NCBI taxonomy database
    filt = df_taxoNCBI['Organism'].isin(organism) #new is a boolean series indicating if organism in taxonomy is in that of df_phylo
    df_tax_sele = df_taxoNCBI[filt][['TaxID','Organism','Tax_kingdom','Tax_group']] #extract those True with collum 'ID','Organism','Tax_kingdom','Tax_group'
    df_tax_sele['Organism'] = df_tax_sele['Organism'].str.strip().str.replace(' ', '_')
    df_seq_tax = pd.merge(df_tax_sele,df_sum, on='Organism')# join the df_tax_sele and df_sum based on the common collumn 'organism'
    #using efetch to get the taxonomy of blastname (more common)
    tax_ID = list(df_seq_tax['TaxID'])
    Tax_blastName = []
    for id in tax_ID:
        do_efetch = 'efetch -db taxonomy -id ' + str(id)
        result = sp.run(do_efetch, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True)
        Tax_blastName.append(result.stdout.decode("utf-8").strip().split('\n')[-1].split(',')[-1])
    df_seq_tax['Tax_blastName'] = pd.Series(Tax_blastName).str.strip().str.replace(' ', '_')
    # df_seq_tax.drop_duplicates(inplace=True)
    #export to a file 'hit_summary_taxo.csv': ID,Organism,Tax_kingdom,Tax_group,Accession,Description,E-value,Identity,Sequence,Tax_blastName
    df_seq_tax.to_csv('hit_summary_taxo.csv', index=False)
    df_seq_tax['Sequence'].to_csv('homologous_sequences.fasta', index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ")
##############################################################################################################################
# taxo()
