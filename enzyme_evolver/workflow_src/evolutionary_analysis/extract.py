#extract single tax (cyanobacteria as an example)
import subprocess as sp
import pandas as pd
import csv
def sequence_process(sequence_tax):
    
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
            #do multiple sequence alignment
            muscle = 'muscle -in ' + outfile + ' -out ' + each +'_muscle.fasta'
            sp.run(muscle, shell=True)
            #mveiw and pdf for initial sequence alginment
            mview = 'mview ' + each +'_muscle.fasta -threshold 100  -consensus on -con_threshold 90  -con_ignore class -con_coloring identity -html head -bold -css on -coloring group  > ' + each + '_mview.html'
            xml2pdf = 'wkhtmltopdf  -s B0  ' + each + '_mview.html ' + each + '_mview.pdf'
            sp.run(mview,shell=True)
            sp.run(xml2pdf,shell=True)
            ##Trim the sequence alignment by using GBlocks
            gblock = 'Gblocks '  + each +'_muscle.fasta' + ' -t=p -e=-gb -b4=10'
            sp.run(gblock, shell=True)
            #color the conserved motif by mview and save as pdf by wkhtmltopdf
            mview_gblock = 'mview ' + each +'_muscle.fasta-gb -threshold 100  -consensus on -con_threshold 90  -con_ignore class -con_coloring identity -html head -bold -css on -coloring group  > ' + each + '_mview_gblock.html'
            xml2pdf_gblock = 'wkhtmltopdf  -s B0  ' + each + '_mview_gblock.html ' + each + '_mview_gblock.pdf'
            sp.run(mview_gblock,shell=True)
            sp.run(xml2pdf_gblock,shell=True)
sequence_tax = 'seq_tax_refine.csv'
sequence_process(sequence_tax)
