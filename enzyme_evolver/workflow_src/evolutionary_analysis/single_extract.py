#Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com
####################################################################################################################################################################
#function: extract sequences within single family (eg.cyanobacteria)    
#Input: hit_summary_taxo.csv from tax.py                                                                                                            
#Usage: python3 single_extract.py
#output file: eg. 0. cyanobacteria_seq.fasta, 1. cyanobacteria_muscle.fasta (sequence alignment by using Muscle); 2. cyanobacteria_mview.pdf (a pdf to view the sequence alignment); 
    # 3. cyanobacteria_muscle.fasta-gb (Trimmed seqeunce alignment file using Gblock ); 4. cyanobacteria_mview_gblock.pdf (a pdf to view Trimmed seqeunce alignment)
    #output:family-seq.fasta eg. cyanobacteria-seq.fasta; sequence alignment file eg. cyanobacteria_muscle.fasta; Trimmed sequence alignment file eg. cyanobacteria_mucle.fasta-gb; and corresponding sequence alignment pdf and images files
#####################################################################################################################################################################

import subprocess as sp
import pandas as pd
import csv
from pdf2image import convert_from_path


def pdf2tif(ifile,ofile):
    pages = convert_from_path(ifile)
    for page in pages:
        page.save(ofile, 'TIFF')


def sequence_process(file):
    #function: 1. sequence alignment; 2. mview 3. Gblock triming; 4. mview the trimed seqeunce alignment file
    #output file: eg: 1. cyanobacteria_muscle.fasta (sequence alignment by using Muscle); 2. cyanobacteria_mview.pdf (a pdf to view the sequence alignment); 
    # 3. cyanobacteria_muscle.fasta-gb (Trimmed seqeunce alignment file using Gblock ); 4. cyanobacteria_mview_gblock.pdf (a pdf to view Trimmed seqeunce alignment)
    #do multiple sequence alignment
    each = file.split('.')[0]
    muscle = 'muscle -in ' + file + ' -out ' + each +'_muscle.fasta'
    sp.run(muscle, shell=True)
    #mveiw and pdf for initial sequence alginment
    #mview cyanobacteria_muscle.fasta -threshold 100  -consensus on -con_threshold 80  -con_ignore class  -con_coloring identity -html head -bold -css on -coloring consensus  -conservation on> mview.html
    mview = 'mview ' + each +'_muscle.fasta -threshold 100  -consensus on -con_threshold 80  -con_ignore class  -con_coloring identity -html head -bold -css on -coloring consensus  -conservation on  > ' + each + '_mview.html'
    xml2pdf = 'wkhtmltopdf  -s B0  ' + each + '_mview.html ' + each + '_mview.pdf'
    sp.run(mview,shell=True)
    sp.run(xml2pdf,shell=True)
    pdf2tif(each + '_mview.pdf',each + '_mview.tif')
    ##Trim the sequence alignment by using GBlocks
    gblock = 'Gblocks '  + each +'_muscle.fasta' + ' -t=p -e=-gb -b4=10'
    sp.run(gblock, shell=True)
    #color the conserved motif by mview and save as pdf by wkhtmltopdf
    mview_gblock = 'mview ' + each +'_muscle.fasta-gb -threshold 100  -consensus on -con_threshold 80  -con_ignore class  -con_coloring identity -html head -bold -css on -coloring consensus  -conservation on  > ' + each + '_mview_gblock.html'
    xml2pdf_gblock = 'wkhtmltopdf  -s B0  ' + each + '_mview_gblock.html ' + each + '_mview_gblock.pdf'
    sp.run(mview_gblock,shell=True)
    sp.run(xml2pdf_gblock,shell=True)
    pdf2tif(each + '_mview_gblock.pdf',each + '_mview_gblock.tif')

def single_align(infile='hit_summary_taxo.csv'):
    #function: extract sequences in one family (using blastname in NCBI) and sequence alignment within the single family
    #output file: eg. 0. cyanobacteria_seq.fasta, 1. cyanobacteria_muscle.fasta (sequence alignment by using Muscle); 2. cyanobacteria_mview.pdf (a pdf to view the sequence alignment); 
    # 3. cyanobacteria_muscle.fasta-gb (Trimmed seqeunce alignment file using Gblock ); 4. cyanobacteria_mview_gblock.pdf (a pdf to view Trimmed seqeunce alignment)
    #output:family-seq.fasta eg. cyanobacteria-seq.fasta; sequence alignment file eg. cyanobacteria_muscle.fasta; Trimmed sequence alignment file eg. cyanobacteria_mucle.fasta-gb; and corresponding sequence alignment pdf and images files

    df_seq_tax = pd.read_csv(infile)
    #turn space to _,eg. red algae as red_algae
    df_seq_tax['Tax_blastName'] = df_seq_tax['Tax_blastName'].str.strip().str.replace(' ', '_')
    #how many families? what are they?
    tax = df_seq_tax['Tax_blastName'].unique()
    #save to a single sequence file for each family
    for each in list(tax):
        #filter the data by tax_blastname and save it a new dataframe df_each
        filt = (df_seq_tax['Tax_blastName'] ==each)
        df_each = df_seq_tax[filt]
        print(df_each)
        #how many sequences in the family
        seq_number = len(df_each)
        # print(seq_number)
        #dumpy seqeunces in single family to output file
        outfile = each+'.fasta'
        df_each['Sequence'].to_csv(outfile, index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ", header=False)
        #sequence processing for each family which containing over 2 sequences
        if seq_number > 1:
            sequence_process(outfile)

#############################################################################################################################################
single_align()