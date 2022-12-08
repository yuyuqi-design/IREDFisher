#Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com
####################################################################################################################################################################
#function: summarize blast details      
#Input: blast.csv from(run_blast.py) protein_name: eg.Halogenase                                                                                                                
#Usage: python3 hit_summary.py blast.csv Halogenase
#output: hit_summary.csv; sequence.fasta (all sequences in blast.csv) which can be used for sequence alignment.
#####################################################################################################################################################################
import sys
import subprocess as sp
import pandas as pd
import csv
from workflow_src.evolutionary_analysis.get_fasta import get_fasta
#requir csv file which is easier to process by pandas
def hit_summary(df=None,file='', protein_name='', save_csv=True):
        if file == '':
            df_sum = df
        else:
            df_sum = pd.read_csv(file)
        #print(df_sum)
        #extract the Accession column and for each item, split it with '.' as delimiter; 
        #with expand = True specify that the split item is sepreated into two columns
        #accession = (df_phylo['Accession']).str.split('.', expand = True)
        #save the first column from splitted columns
        accesion_number = pd.Series(df_sum['Accession'])
        #write accesion_number into a file as csv without index(index = False)
        accesion_number.to_csv('accesion_number.txt', index = False,header=False)
        #download fasta file from ncbi
        get_fasta('accesion_number.txt',ofolder)

        #make a pandas series based on sequence data in accesion number
        #all_seq is a list with each element for each sequence
        all_seq = []
        for code in list(accesion_number):
                with open(code + '.fasta') as f:
                        #read the content
                        all = f.readlines()
                        #extract the first line
                        firs_line = all[0]
                        #extract the organism from first line
                        orsm= '_'.join(firs_line.strip().split(']')[-2].split('[')[-1].split())
                        #extract the sequence
                        seq = all[1:]
                        #join each line in fasta file into one line to make a one-line sequece
                        seq = ''.join(all[1:])
                        #combine the organism, description and sequence
                #         if ('reductase A' in firs_line) or ('reductase 1' in firs_line):
                #             orsm_des = orsm + '_PORA'
                #         elif ('reductase B' in firs_line) or ('reductase 2' in firs_line):
                #             orsm_des = orsm + '_PORB'
                #         elif ('reductase C' in firs_line) or ('reductase 3' in firs_line):
                #             orsm_des = orsm + '_PORC'
                #         elif ('isoform X1' in firs_line):
                #             orsm_des = orsm + '_POR_isoform_X1'
                #         elif ('isoform X2' in firs_line):
                #             orsm_des = orsm + '_POR_isoform_X2'
                #         else:
                #             orsm_des = orsm + '_POR'
                        orsm_des = orsm + '' + protein_name.strip()
                        orsm_seq = '>' +orsm_des + '\n' + seq 
                #         orsm_seq = '>' +orsm + '\n' + seq
                        #append each acceaccesion_number sequence to all_seq
                        all_seq.append(orsm_seq)
        # description column: protochlorophyllide reductase [Thermosynechococcus elongatus]
        #split description column into two seperated columns: description for enzyme name, organism for the cell type
        # description= pd.Series(df_phylo['Description']).str.split('[', expand = True)[0]

        # organism = pd.Series(df_phylo['Description']).str.split('[', expand = True)[1].str.split(']', n = 1, expand = True)[0]

        # #covert all_seq (list) to pandas series
        # sequence = pd.Series(all_seq)
        # #use pandas concat to join different series vertically(column1,column2,...) by specifying axis = 1
        # all_join = pd.concat([accesion_number,description, organism,sequence], axis = 1)
        #write joined series to csv file, with header as ['Accession', 'Description', 'Organism', 'Sequence']
        #Accession,Description, Organism,E, Identity
        df_sum['Sequence'] = all_seq
        # delete all the fasta files and clean the directory
        #sp.run('rm *fasta', shell = True)

        if save_csv == True:
                df_sum.to_csv('hit_summary.csv', index=False)
        else:
                return df_sum


#########################################################################################################################################
# file=sys.argv[1]
# protein_name = sys.argv[2]
# hit_summary(file,protein_name)
