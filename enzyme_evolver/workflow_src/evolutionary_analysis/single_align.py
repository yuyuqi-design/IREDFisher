#extract single tax (cyanobacteria as an example)
from sequence_process import sequence_process
import subprocess as sp
import pandas as pd
import csv
def single_align(sequence_tax):
    #function: extract sequences in one family (using blastname in NCBI) and sequence alignment within the single family
    #output file: eg. 0. cyanobacteria_seq.fasta, 1. cyanobacteria_muscle.fasta (sequence alignment by using Muscle); 2. cyanobacteria_mview.pdf (a pdf to view the sequence alignment); 
    # 3. cyanobacteria_muscle.fasta-gb (Trimmed seqeunce alignment file using Gblock ); 4. cyanobacteria_mview_gblock.pdf (a pdf to view Trimmed seqeunce alignment)

    df_seq_tax = pd.read_csv(sequence_tax)
    #turn space to _,eg. red algae as red_algae
    df_seq_tax['Tax_blastName'] = df_seq_tax['Tax_blastName'].str.strip().str.replace(' ', '_')
    #how many families?
    tax = df_seq_tax['Tax_blastName'].unique()
    #save to a single sequence file for each family
    for each in list(tax):
        #filter the data by tax_blastname and save it a new dataframe df_each
        filt = (df_seq_tax['Tax_blastName'] ==each)
        df_each = df_seq_tax[filt]
        #how many sequences in the family
        seq_number = len(df_each)
        #dumpy seqeunces in single family to output file
        outfile = each+'_seq.fasta'
        df_each['Sequence'].to_csv(outfile, index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ", header=False)
        #sequence processing for each family which containing over 2 sequences
        if seq_number > 1:
            sequence_process(outfile))