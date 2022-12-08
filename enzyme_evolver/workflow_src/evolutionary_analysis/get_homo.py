import subprocess as sp
import uuid
import pandas as pd
import sys
import csv
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from pathlib import Path
from enzyme_evolver.workflow_src.evolutionary_analysis import get_fasta
# require cd-hit installed
def run_blast_remote(query, refseqdb='nr', ofolder=None,max_seq='400'):
    # Blast in NCBI refseq_protein database or in-house database(upload by users)
    # input: query
    # output: a xml file
    # -outfmt equals 5 means an xml format result
    # option1: specify the protein name (eg, halogenase)
    # option2 : one sequence one organism
    # option3: upload an in-house database (later) by converting fasta format to database format using
    #         makeblastdb -in inhouse.fasta -parse_seqids -blastdb_version 5 -title "inhouse" -dbtype prot
    #############################################blast remotely
    # fasta_string = open(query).read()
    # result_handle = NCBIWWW.qblast('blastp', refseqdb, fasta_string,hitlist_size=400, results_file='blast.xml')
    command = 'blastp' + ' -query ' + query + ' -db ' + refseqdb + ' -out ' + ofolder+ 'blast.xml'+ ' -outfmt 5 -max_hsps 1 -max_target_seqs '+ max_seq+ ' -remote'
    sp.run(command, shell=True)

def parse_blast(protein='', file='blast.xml', save_csv=True, ofolder=''):
    #Parse the blast result
    #input: blast xml result; the name of the enzyme. for example: halogenase
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
            if len(alignment.title.strip().split('[')) < 2:
                pass
            else:
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
                    #print(hsp.num_alignments)
                    identity = hsp.identities/hsp.align_length
                    #if (identity > 0.3) and (hsp.expect < 0.01):
                    E.append(hsp.expect)
                    Identity.append(identity)
    #print(len(E))
    #print(len(Organism))
    result_handle.close()
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
    df_unqOrgm = df_unqOrgm.dropna()

    #Header: Accession,Description, Organism,E, Identity
    if save_csv == True:
        df_unqOrgm.to_csv(ofolder+'blast.csv',index=False)
    else:
        return df_unqOrgm

def hit_summary(df=None, file='', protein_name='', save_csv=True, ofolder=None):
    if file == '':
        df_sum = df
    else:
        df_sum = pd.read_csv(file)
    # print(df_sum)
    # extract the Accession column and for each item, split it with '.' as delimiter;
    # with expand = True specify that the split item is sepreated into two columns
    # accession = (df_phylo['Accession']).str.split('.', expand = True)
    # save the first column from splitted columns
    accesion_number = pd.Series(df_sum['Accession'])
    # write accesion_number into a file as csv without index(index = False)
    accesion_number.to_csv(ofolder+'accesion_number.txt', index=False, header=False)
    # download fasta file from ncbi
    get_fasta.get_fasta(ofolder+'accesion_number.txt',ofolder)

    # make a pandas series based on sequence data in accesion number
    # all_seq is a list with each element for each sequence
    all_seq = []
    for code in list(accesion_number):
        with open(ofolder+code + '.fasta') as f:
            #print(ofolder+code + '.fasta')
            # read the content
            all = f.readlines()
            #print(f.readlines())
            #print(all)
            # extract the first line
            firs_line = all[0]
            # extract the organism from first line
            orsm = '_'.join(firs_line.strip().split(']')[-2].split('[')[-1].split())
            # extract the sequence
            seq = all[1:]
            # join each line in fasta file into one line to make a one-line sequece
            seq = ''.join(all[1:])
            # combine the organism, description and sequence
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
            orsm_seq = '>' + orsm_des + '\n' + seq
            #         orsm_seq = '>' +orsm + '\n' + seq
            # append each acceaccesion_number sequence to all_seq
            all_seq.append(orsm_seq)
    # description column: protochlorophyllide reductase [Thermosynechococcus elongatus]
    # split description column into two seperated columns: description for enzyme name, organism for the cell type
    # description= pd.Series(df_phylo['Description']).str.split('[', expand = True)[0]

    # organism = pd.Series(df_phylo['Description']).str.split('[', expand = True)[1].str.split(']', n = 1, expand = True)[0]

    # #covert all_seq (list) to pandas series
    # sequence = pd.Series(all_seq)
    # #use pandas concat to join different series vertically(column1,column2,...) by specifying axis = 1
    # all_join = pd.concat([accesion_number,description, organism,sequence], axis = 1)
    # write joined series to csv file, with header as ['Accession', 'Description', 'Organism', 'Sequence']
    # Accession,Description, Organism,E, Identity
    df_sum['Sequence'] = all_seq
    # delete all the fasta files and clean the directory
    # sp.run('rm *fasta', shell = True)
    get_fasta.clean_fasta(ofolder+'accesion_number.txt',ofolder)
    df_sum['Sequence'].to_csv(ofolder+'homologous_sequences.fas', index=False, quoting=csv.QUOTE_NONE, quotechar="", escapechar=" ",
                        header=False)
    if save_csv == True:
        df_sum.to_csv(ofolder+'hit_summary.csv', index=False)

    return df_sum #['Accession', 'Description', 'Organism', 'Sequence']

###cluster the fasta sequences by threshold 90%

def seqcluster(seq='homologous_sequences.fas',ofolder=None):
    seq = str(seq).strip()
    print(seq.split('.')[0])
    command_cluster = 'cd-hit -i ' + ofolder+ seq + ' -o ' + ofolder+ seq.split('.')[0] + ' -T 2 -c 0.75'
    #print(command_cluster)
    sp.run(command_cluster, shell=True)
    # command_getRep = "grep '\*' " + seq.split('.')[0]+".clstr|awk -F \> '{print $2}'|awk -F . '{print $1}'|sort > seq_rep_id_" + seq.split('.')[0]+".txt"
    # sp.run(command_getRep, shell=True)
    # command_catsequence = "sh catsequence90.sh " + "seq_rep_id_" + seq.split('.')[0]+".txt"
    # sp.run(command_catsequence, shell=True)
    command_cluster_organism = "grep \> " + ofolder+"homologous_sequences |awk -F \> '{print $2}' > " +  ofolder+"cluster_organism.txt"
    sp.run(command_cluster_organism, shell=True)
    #sp.run('rm *fasta', shell=True)

#function: remove the low quality sequences
#Input: hit_summary_taxo.csv       dump_name.txt
#Usage: python3 sequence_dump.py
#output: hit_summary_taxo_refine.csv
#####################################################################################################################################################################
def sequence_dump(df=None, file='', seq_dump="cluster_organism.txt"):

    if df is None:
        df = pd.read_csv(file)

    #seq_tax is a file containing all information: ID,Organism,Tax_kingdom,Tax_group,Unnamed: 0,Accession,Description,Sequence,Tax_blastName;
    #seq_dump contains the organism name of low quality sequences eg:Helianthus_annuus_PORB
    with open(seq_dump) as file:
        dump = file.readlines()
        dump = [i.strip() for i in dump]
        #print(dump)
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
        #df_refine.to_csv('hit_summary.csv' + fileid, index=False)
        #df_refine['Sequence'].to_csv('homologous_clustered90.fas' + fileid, index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ", header=False)

        return df_refine

# This function was used to generate taxoNCBI.csv
def load_dftaxoNCBI(path='enzyme_evolver/database/fullnamelineage.dmp'):
    print('LOADING FULLLINEAGE.DMP...')
    df_taxoNCBI = pd.read_csv(path, delimiter='|', names = ['TaxID','Organism','Tax','nan'])
    print('FULLLINEAGE.DMP LOADED OK')

    # remove tab from organism and tax collum
    print('PARSING DF..')
    df_taxoNCBI['Organism'] = df_taxoNCBI['Organism'].str.strip()
    df_taxoNCBI['Tax_kingdom'] = df_taxoNCBI['Tax'].str.strip().str.split(';', expand=True)[1]
    df_taxoNCBI['Tax_group'] = df_taxoNCBI['Tax'].str.strip().str.split(';', expand=True)[4]
    print('DF PARSED')
    df_taxoNCBI.to_csv('taxoNCBI.csv')


def taxo(df=None,file='',
         save_csv=True,
         ofsummary='hit_summary.csv',
         ofseq='homologous_sequences90.fas',
         taxoNCBI='enzyme_evolver/database/taxoNCBI.csv'):
    if file == '':
        df_sum = df
    else:
        #get the organism name from selected seqeunces
        df_sum = pd.read_csv(file)

    organism = df_sum['Organism'].str.strip().str.replace('_', ' ')
    #print(organism)
    #NCBI taxonomy database: taxID, organism,taxonomy, nan

    df_taxoNCBI = pd.read_csv(taxoNCBI)

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
    df_seq_tax.to_csv(ofsummary, index=False)
    df_seq_tax['Sequence'].to_csv(ofseq, header=False, index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ")
    return df_seq_tax

