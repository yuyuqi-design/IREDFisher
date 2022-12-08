#Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com
####################################################################################################################################################################
#function: sequence alignment for all sequences       
#Input: hit_summary_taxo_refine.csv                                                                                                                  
#Usage: python3 msa.py
#output: hit_refine.fasta; hit_refine_muscle.fasta,hit_refine_muscle.fasta-gb
#####################################################################################################################################################################
import sequence_process
import pandas as pd
import csv
#df = pd.read_csv('hit_summary_taxo_refine.csv')
#df.sort_values(by='Tax_blastName',ascending=True, inplace=True)
#df['Sequence'].to_csv('hit_refine.fasta', index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ", header=False)
sequence_process('homologous_clustered90.fas')