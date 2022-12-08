#Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com
####################################################################################################################################################################
#function: remove the low quality sequences       
#Input: hit_summary_taxo.csv       dump_name.txt                                                                                                         
#Usage: python3 sequence_dump.py
#output: hit_summary_taxo_refine.csv
#####################################################################################################################################################################
import pandas as pd
import csv
def sequence_dump(seq_tax='hit_summary_taxo.csv', seq_dump='cluster_organism.txt'): 

    #seq_tax is a file containing all information: ID,Organism,Tax_kingdom,Tax_group,Unnamed: 0,Accession,Description,Sequence,Tax_blastName;
    #seq_dump contains the organism name of low quality sequences eg:Helianthus_annuus_PORB
    df = pd.read_csv(seq_tax)
    with open(seq_dump) as file:
        dump = file.readlines()
        dump = [i.strip() for i in dump]
        print(dump)
        #print(dump)
        # dump = '|'.join(dump)
#        print(dump)
#         print(dump)
        #filt = (df['Organism'] == dump)
        df_uniq = df.drop_duplicates(subset='Organism', keep='first')
        df_refine = df_uniq[df_uniq['Organism'].isin(dump)]
        # print(len(df))
        # print(len(df_refine))
        #filt = df['Organism'] == dump
        #print(df['Organism'])
        df_refine.to_csv('hit_summary_taxo_refine.csv', index=False)
        df_refine['Sequence'].to_csv('homologous_clustered90.fas', index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ", header=False)
        #return df[filt]
#################################################################################################
# sequence_dump()
# seq_tax = 'seq_tax.csv'
# seq_dump = 'seq_dump.dat'
# sequence_dump(seq_tax,organism_file)
    