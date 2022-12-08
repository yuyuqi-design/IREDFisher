# coding: utf-8
# Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com

#######################################################################################################################################
#                                                                                                                                     #
#                                   This program is called Auto_modeller.                                                             # 
#                                                                                                                                     #
# It is constructed to help experimental researchers get three dimentional structures with higher working effieciency. Only by giving #
# a Uniprot/NCBI code or a series of code saved in a text file, you will automatically get the sequence of your protein in FASTA      #
# format. Then the program will check in PDB non-redundant database if there is already a experimental structure. If yes, the program #
# will give you the experimental structure in a new directory with the name Experi_Stru_$ID.                                          #
# Otherwise, Auto_modeller will make sequence alignment and find a template to build a model of your protein by using the program of  #
# Modeller 9.23.                                                                                                                      #
# Note: only sequence identity between user's sequence and template protein's sequence over 30% will be build.                        #
# Usage: ./auto_modeller.py file.txt > log  or ./auto_modeller.py P12345 > log                                                         #
#                                                                                                                                     #
#######################################################################################################################################

import os
from urllib import request
import sys
import subprocess as sp
from modeller import *
import pandas as pd
from modeller.automodel import *
from modeller.scripts import complete_pdb
from Bio import Entrez,SeqIO, AlignIO
import shutil


def fasta2pir(seqcode):
    ###############################################################################
    #This function is to convert the sequence from FASTA to PIR format which is   #
    #recognized by Modeller                                                       #
    ###############################################################################'
    e = environ()
    #convert fasta to pir format
    a = alignment(e, file=seqcode + '.fasta', alignment_format='FASTA')
    a.write(file=str(seqcode).strip() + '_0.pir', alignment_format='PIR')
    #remove the first empty line
    tempofile = str(seqcode).strip() + '_0.pir'
    filename = str(seqcode).strip() + '.pir'
    with open(tempofile) as f:
        with open(filename, 'w') as f2:
            seq = f.readlines()[3:]
            #print(seq)
            f2.write('>P1;' + str(seqcode).strip()+'\n')
            line = 'sequence:' + str(seqcode).strip() + ':::::::0.00: 0.00'
            f2.write(line + '\n')
            f2.writelines(seq)
            #remov = "sed -i '/^$/d' " + str(code).strip() + '.pir'
    #sp.run(remov, shell=True)

def template_search(code,db='enzyme_evolver/database/pdb_95.pir'):
    ###################################################################################
    #This function is to make sequence alignment between user's sequence and proteins #
    # in PDB database (pdb_95.pir which is updated in the end of November, 2019)      #
    ###################################################################################

    #    global template, template_chain, template_PdbCode, Experi_PdbCode

    #log.verbose()
    env = environ()

    #-- Prepare the input files

    #-- Read in the sequence database
    sdb = sequence_db(env)
    sdb.read(seq_database_file=db, seq_database_format='PIR',
     chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

    #-- Write the sequence database in binary form
    sdb.write(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
      chains_list='ALL')

    #-- Now, read in the binary database
    sdb.read(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
     chains_list='ALL')

    #-- Read in the target sequence/alignment
    aln = alignment(env)
    aln.append(file=str(code).strip() + '.pir', alignment_format='PIR', align_codes='ALL')

    #-- Convert the input sequence/alignment into
    #   profile format
    prf = aln.to_profile()

    #-- Scan sequence database to pick up homologous sequences
    prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
      gap_penalties_1d=(-500, -50), n_prof_iterations=1,
      check_profile=False, max_aln_evalue=0.01)

    #-- Write out the profile in text format
    prf.write(file=str(code).strip() + '.prf', profile_format='TEXT')

    #-- Convert the profile back to alignment format
    aln = prf.to_alignment()

    #-- Write out the alignment file
    aln.write(file=str(code).strip() + '.ali', alignment_format='PIR')

    #######################################################################################

def template_des(code):
    code=code.strip()
    with open(code) as f:
        # read the content
        all = f.readlines()
        # extract the first line
        firs_line = all[0]
        # extract the description from first line
        seq_details = firs_line.strip().split('|')
        organism = firs_line.strip().split('|')[3]
        description = firs_line.strip().split('|')[2]
    return description, organism



def check_template(code):
    ###################################################################################
    #This function is to select best template to build homology model                 #
    ###################################################################################
    # global template, template_chain, template_PdbCode, Experi_PdbCode
    # global templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain, template2_PdbCode, template3_chain, template3_PdbCode
    cwd = os.getcwd()
    #prf including the alignment, sequence identiy, E-vaule etc information
    result_file = str(code).strip() + '.prf'
    #read potential template file from prf file
    df = pd.read_csv(result_file, skiprows=6, delim_whitespace = True, names=["index", "PDB", "SX", "start_template", 'end_template','seq_startmatch', 'seq_endmatch', 'template_startmatch','template_endmatch','Cover','SI','E','ALI'])

    #filter the sequence identity =100 and E-value = 0 so can check if there is a experimental structure in structure database.
    Experi_PdbCode = list(df[(df['SI'] == 100.0) & (df['E'] == 0.0)]['PDB'])
    if Experi_PdbCode: #if Experi_PdbCod is not empty
    #make a directory and save the experimental structure into that directory
    #sp.run('mkdir ' + 'Experi_Stru_'+ str(code).strip(), shell = True)
        os.mkdir('Experi_Stru_'+ str(code).strip())
        path = cwd + '/Experi_Stru_'+ str(code).strip()
        #os.chdir(path)
        #sp.run('cd ' + 'Experi_Stru_'+ str(code).strip(), shell = True)
        for i in Experi_PdbCode:
        #get_pdb = 'wget http://www.rcsb.org/pdb/files/' + str(i[0:4]).strip() + '.pdb'
        #sp.run(get_pdb, shell = True)
            request.urlretrieve('http://files.rcsb.org/download/' + str(i[0:4]).strip() +'.pdb', 'Experi_Stru_'+ str(code).strip()+'/'+ str(i[0:4]).strip() +'.pdb')
            print('This sequence has experimental structure. Congratulations! :' + str(i[0:4]) + '.pdb' + '\n' + 'DONE!')
            try:
                shutil.copy(str(i[0:4]).strip() +'.pdb', cwd)
            except:
                pass
            #os.chdir(cwd)
    #else:
    df_sort = df.sort_values(by='SI', ascending=False)
    #check if the biggest sequence identiry is over 30%
    biggest_SI = float(df_sort.iloc[0,10])
    if biggest_SI > 0.0:
    #example is template:2aqjA; template_PdbCode:2aqj; template_chain : A
    #single template modelling
        # templates10 = list(df_sort.iloc[:10, 1])
        # df_temp = pd.DataFrame(templates10)
        # df_temp['code'] = df_temp.iloc[:, 0].str[0:4]
        # df_temp['chain'] = df_temp.iloc[:, 0].str[-1]
        # df_temp.drop_duplicates(subset=['code'], inplace=True)
        # templates = list(df_temp.iloc[:3, 0])

        df_code_seqID = df_sort[["PDB", "SI",'seq_startmatch', 'seq_endmatch', 'template_startmatch','template_endmatch','Cover']][:10]
        # print(df_code_seqID)
        df_code_seqID['code'] = df_code_seqID['PDB'].str[0:4]
        df_code_seqID['chain'] = df_code_seqID['PDB'].str[-1]
        df_code_seqID.drop_duplicates(subset=['code'], inplace=True)
        df_info = df_code_seqID[:1][['PDB','SI','seq_startmatch', 'seq_endmatch', 'template_startmatch','template_endmatch','Cover']]
        df_info.insert(0,'SEQ',code.strip())
        df_info['Model'] = code.strip()+'.pdb'
        templates = list(df_code_seqID.iloc[:3, 0])
        template1 = str(df_sort.iloc[0, 1])
        template1_PdbCode = template1[0:4]
        template1_chain = template1[-1]
        if len(templates)>=2:
            template2 = templates[1]
            template2_PdbCode = template2[0:4]
            template2_chain = template2[-1]
        else:
            template2 = ''
            template2_PdbCode = ''
            template2_chain = ''
        if len(templates) >= 3:
            template3 = templates[2]
            template3_PdbCode = template3[0:4]
            template3_chain = template3[-1]
        else:
            template3 = ''
            template3_PdbCode = ''
            template3_chain = ''
        print('The template structure for ' + str(code).strip() + 'is:' + str(template1) + '\n')
        download_fasta = f'wget https://www.rcsb.org/fasta/entry/' + template1_PdbCode
        sp.run(download_fasta, shell=True)
        description, organism = template_des(template1_PdbCode)
        df_info['description'] = description
        df_info['organism'] = organism
        df_info.to_csv('info.txt', index=False, header=['Sequence',"template_PDBcode", "SI",'seq_startmatch', 'seq_endmatch', 'template_startmatch','template_endmatch','Cover','Model', 'description','organism'], sep='\t')
        sp.run(f'rm {template1_PdbCode}',shell=True)

        try:
            request.urlretrieve('http://files.rcsb.org/download/' + template1_PdbCode +'.pdb', template1_PdbCode +'.pdb')
            request.urlretrieve('http://files.rcsb.org/download/' + template2_PdbCode +'.pdb', template2_PdbCode +'.pdb')
            request.urlretrieve('http://files.rcsb.org/download/' + template3_PdbCode +'.pdb', template3_PdbCode +'.pdb')
        except:
            pass

    #if no available template
    else:
        templates = ''
        template1 = ''
        template1_PdbCode = ''
        template1_chain = ''
        template2 = ''
        template2_PdbCode = ''
        template2_chain = ''
        template3 = ''
        template3_PdbCode = ''
        template3_chain = ''
        f=open('info.txt', 'w')
        f.write(f'{code}:no reliable template was found\n')
        f.close()
        print('Sorry, there is no reliable template  in the database\n')

    return templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain, template2_PdbCode, template3_chain, template3_PdbCode
    




###########################################################################################
def salign(seqcode,template1_PdbCode='',template1_chain='',template2_PdbCode='', template2_chain='',template3_PdbCode='', template3_chain=''):
    # Illustrates the SALIGN multiple structure/sequence alignment

    env = environ()
    env.io.atom_files_directory = './:../atom_files/'

    aln = alignment(env)
    for (Code, chain) in ((template1_PdbCode, template1_chain ), (template2_PdbCode, template2_chain), (template3_PdbCode, template3_chain)):
        mdl = model(env, file=Code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
        aln.append_model(mdl, atom_files=Code, align_codes=Code+chain)

    for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                        ((1., 0.5, 1., 1., 1., 0.), False, True),
                                        ((1., 1., 1., 1., 1., 0.), True, False)):
        aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                rr_file='$(LIB)/as1.sim.mat', overhang=30,
                gap_penalties_1d=(-450, -50),
                gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                dendrogram_file='dedro.tree',
                alignment_type='tree', # If 'progresive', the tree is not
                                        # computed and all structues will be
                                        # aligned sequentially to the first
                feature_weights=weights, # For a multiple sequence alignment only
                                            # the first feature needs to be non-zero
                improve_alignment=True, fit=True, write_fit=write_fit,
                write_whole_pdb=whole, output='ALIGNMENT QUALITY')

    aln.write(file=str(seqcode).strip() + 'salign.pap', alignment_format='PAP')
    aln.write(file=str(seqcode).strip() + 'salign.ali', alignment_format='PIR')

    aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
            rr_file='$(LIB)/as1.sim.mat', overhang=30,
            gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
            gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
            alignment_type='progressive', feature_weights=[0]*6,
            improve_alignment=False, fit=False, write_fit=True,
            write_whole_pdb=False, output='QUALITY')

######################################################################################


def align2d(code,template1,template1_PdbCode,template1_chain):


   ########################################################################################
   # This function is to align the user's sequence to selected template sequence          #
   ########################################################################################

    env = environ()
    aln = alignment(env)
    mdl = model(env, file=template1_PdbCode, model_segment=('FIRST:'+ template1_chain,'LAST:' + template1_chain))
    #use the template variable in template_search module.
    aln.append_model(mdl, align_codes=template1, atom_files=template1_PdbCode + '.pdb')
    #use the code variable from fasta2pir function
    aln.append(file=str(code).strip() + '.pir', align_codes= str(code).strip())
    aln.align2d()
    aln.write(file=str(code).strip() + '-' + template1 + '.ali', alignment_format='PIR')
    aln.write(file=str(code).strip() + '-' + template1 + '.pap', alignment_format='PAP')

#########################################################################################
def align2d_mult(code):
    env = environ()

    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read aligned structure(s):
    aln = alignment(env)
    aln.append(file=str(code).strip() +  'salign' + '.ali', align_codes='all')
    aln_block = len(aln)

    # Read aligned sequence(s):
    aln.append(file=str(code).strip() + '.pir', align_codes= str(code).strip())

    # Structure sensitive variable gap penalty sequence-sequence alignment:
    aln.salign(output='', max_gap_length=20,
            gap_function=True,   # to use structure-dependent gap penalty
            alignment_type='PAIRWISE', align_block=aln_block,
            feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
            gap_penalties_1d=(-450, 0),
            gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
            similarity_flag=True)

    aln.write(file=str(code).strip() + '-' + 'multi' + '.ali', alignment_format='PIR')
    aln.write(file=str(code).strip() + '-' + 'multi' + '.pap', alignment_format='PAP')
    ########################################################################################

def trim_TerGap(file1,file2='temp.fasta', format1='pir',format2='fasta'):

    aln = SeqIO.parse(file1,format1)
    SeqIO.convert(file1,format1,file2,format2)

    aln = AlignIO.read(file2, format2)
    for col in range(aln.get_alignment_length()):
        if not "-" in aln[:,col]:
            position = col
            break
    for col2 in reversed(range(aln.get_alignment_length())):
        if not "-" in aln[:,col2]:
            position2 = col2+1
            break
    SeqIO.write(aln[:,position:position2],file1,format1)
    sp.run(f"sed -i 's/>XX/>P1/g' {file1}", shell=True)
    #remove the starting residue number and ending residue number in template structure
    with open(file1) as f:
        lines = f.readlines()
        line1 = lines[1].split(':')
        #remove starting residue number
        line1[2] = ' '
        # remove ending residue number
        line1[4] = ' '
        lines[1] = ':'.join(line1)
        with open(file1,'w') as fo:
            fo.writelines(lines)


def build_model(code,template1):
   ########################################################################################
   # This function is to build a model based on the selected templates                    #
   ########################################################################################
    #from modeller import soap_protein_od

    env = environ()
    a = automodel(env, alnfile=str(code).strip() + '-' + template1 + '.ali',
                  knowns=template1, sequence=str(code).strip(),
                  assess_methods=(assess.DOPE,
                                  #soap_protein_od.Scorer(),
                                  assess.GA341))
    a.starting_model = 1
    a.ending_model = 1
    try:
        a.make()
    except:
        pass
    
###########################################################################################
def build_multi_models(code,template1,template2,template3):

    env = environ()
    a = automodel(env, alnfile=str(code).strip() + '-' + 'multi' + '.ali',
                knowns=(template1, template2, template3), sequence=str(code).strip(), assess_methods=(assess.DOPE, assess.GA341))
    a.starting_model = 1
    a.ending_model = 1
    a.make()
##########################################################################################


def single2dimer_ali(code, start_residue,end_residue, N_chains='2'):
    code = code.strip()
    template1 = 'singleT' + 'A'
    template1_PdbCode = 'singleT'
    template1_chain = 'A'
    align2d(code, template1,template1_PdbCode,template1_chain)
    with open(code.strip() + '-singleTA.ali' ) as single:
        with open('multi.ali', 'w') as multi:
            ali = single.readlines()
            for line in ali:
                if line.startswith('structure'):
                    seq_struc_BEGIN = ali.index(line) + 1
                if line.startswith('sequence'):
                    seq_target_BEGIN = ali.index(line) + 1 
                    seq_struc_END = seq_target_BEGIN -2
                    seq_target_END = len(ali)
            structure_sequence = ''.join(ali[seq_struc_BEGIN:seq_struc_END]).replace("*","").strip()
            target_sequence = ''.join(ali[seq_target_BEGIN: seq_target_END]).replace("*","").strip()       
            multi.writelines(ali[0])
            multi.write(ali[1].replace('singleTA','dimerT'))
            multi.write('structureX:dimerT:   ' + start_residue + ':A:' + end_residue + ':B:::-1.00:-1.00\n')
            #multi.write(ali[2].replace("A:::-1.00:-1.00","B:::-1.00:-1.00").replace("singleT.pdb","dimerT"))
            multi.write((structure_sequence + '/' + '\n') * (int(N_chains)-1))
            multi.write(structure_sequence + '*')
            multi.write('\n')
            multi.writelines(ali[seq_target_BEGIN-2:seq_target_BEGIN ])
            multi.write((target_sequence + '/' + '\n') * (int(N_chains)-1))
            multi.write(target_sequence + '*')
    command = f"mview -in pir {code}-singleTA.ali|awk '{{if( $1 == 2) print $7,$8,$9}}' |tail -1"
    #print(command)
    perc = sp.run(command, shell=True, stdout=sp.PIPE).stdout.strip().decode("utf-8")
    text_file = open("info.txt", "w")
    text_file.write(perc + '\n')
    text_file.close()
# Demonstrates how to build multi-chain models, and symmetry restraints
#log.verbose()

# Override the 'special_restraints' and 'user_after_single_model' methods:
class MyModel(automodel):
    def special_restraints(self, aln):
        # Constrain the A and B chains to be identical (but only restrain
        # the C-alpha atoms, to reduce the number of interatomic distances
        # that need to be calculated):
        s1 = selection(self.chains['A']).only_atom_types('CA')
        s2 = selection(self.chains['B']).only_atom_types('CA')
        self.restraints.symmetry.append(symmetry(s1, s2, 1.0))
    def user_after_single_model(self):
        # Report on symmetry violations greater than 1A after building
        # each model:
        self.restraints.symmetry.report(1.0)
def build_chains(code):

    env = environ()
    # directories for input atom files
    env.io.atom_files_directory = ['.', '../atom_files']

    # Be sure to use 'MyModel' rather than 'automodel' here!
    a = MyModel(env,
                alnfile  = 'multi.ali' ,     # alignment filename
                knowns   = 'dimerT',              # codes of the templates
                sequence = str(code).strip())              # code of the target
    
    a.starting_model= 1                # index of the first model
    a.ending_model  = 1                # index of the last model
                                    # (determines how many models to calculate)
    a.make()                           # do comparative modeling

######################################################################################
def procheck(code):
    code = code.strip()
    evaluate_cmd = f"procheck {code}.pdb 2.0"
    sp.run(evaluate_cmd,shell=True)
    ##remove files
    rm_cmd = f"rm procheck.prm *.ps *.sco *lot.log *out *sdh *.lan *.pln *sco nb.log *.nb *.rin secstr.log *new anglen.log clean.log"
    sp.run(rm_cmd,shell=True)

def rama_core(code):
    with open(f'{code}.sum') as fin:
        result = fin.readlines()
        for line in result:
            if '| Ramachandran plot:' in line:
                core = float(line.split()[3].strip('%'))
                print(core)


       

