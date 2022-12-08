#requirement:eBoxSize.pl
#Usage: python3 autodocking.py receptors_file.txt ligands_file.txt
#receptors_file writing: 
#4agd.pdb
#4agc.pdb
#...
#ligands_file writing:
#lig1.pdb
#lig2.pdb
#...
#prepare_receptor4.py -r 4agd.pdb -A bonds_hydrogens -e
#prepare_ligand4.py -l lig_4agd.pdb

import os
import subprocess as sp

import numpy as np
from pymol import cmd

from enzyme_evolver.workflow_src.binding_mode.docking import find_site


def mk_dir(rec_file,lig_file,folder):
    #make the sub directory
    name_dir = str(rec_file).split('.')[0] + '_' + str(lig_file).split('.')[0]
    try:
        os.mkdir(folder+name_dir)
    #name the input configure file as rec_lig_config.txt; for example: 4agd_b49_config.txt
    except FileExistsError:
        print('directory exists.')
    return name_dir
#####################################################################################
def align(recs,folder):
    #align all the proteins in receptors_file

    refer_pro = recs[0].strip().split('.')[0]
    print(recs)
    for protein in recs:
        print(folder+protein.strip())
        cmd.load(folder+protein.strip())
    cmd.alignto(refer_pro)
    for protein in recs:
        name = protein.strip().split('.')[0]
        cmd.save(folder+'aln_' + name + '.pdb', name)
    cmd.delete('all')
#####################################################################################
#####################################################################################
def prepare(rec_file, lig_file,folder):
#     global all_ligs
#     for protein in all_proteins:
        #prepare receptors
    # print(os.getcwd())
    #python path: /home/g02808yy/software/mgltools/MGLTools-1.5.6/MGLToolsPckgs/:${PYTHONPATH}
    comand_rec_prepare = 'python2.5 enzyme_evolver/workflow_src/binding_mode/docking/prepare_receptor4.py' \
                         + ' -r ' + folder + 'aln_' + rec_file.strip().split('.')[0]+'.pdb' + '    -A checkhydrogens -U deleteAltB -o ' + folder +str(rec_file).split('.')[0] +'.pdbqt'
    print(comand_rec_prepare)

        #prepare ligands
#     with open(ligands_file) as ligs_file:
#         all_ligs = ligs_file.readlines()
#         for lig in all_ligs:
    comand_lig_prepare = 'python2.5 enzyme_evolver/workflow_src/binding_mode/docking/prepare_ligand4.py -l ' + folder + str(lig_file).split('.')[0]+'.mol2' + ' -o '\
                         + folder+str(lig_file).split('.')[0] +'.pdbqt'
    #print(os.getcwd())
    print(comand_lig_prepare)


    try:
        os.system(comand_rec_prepare)
        pdbqt = folder +str(rec_file).split('.')[0] +'.pdbqt'
        sp.run(f'cp {pdbqt} /home/g02808yy/data/webserver/IREDFisher/enzyme_evolver/database/Public', shell=True)
        os.system(comand_lig_prepare)
    except:
        pass


#####################################################################################
def box_center_residues(rec_file, lig_file, residue_ID1, residue_ID2, residue_ID3, residue_ID4):
    global box_center
    
    with open(rec_file) as receptor: 
        content = receptor.readlines()
        for line in content:
            if (line.strip()[22:27].strip() == str(residue_ID1).strip())  & (line.strip()[13:15] == 'CA') & (line.strip().startswith('ATOM')):
                cord1 = line.strip()[32:55].split()
#                     print(cord1)
            elif (line.strip()[22:27].strip() == str(residue_ID2).strip()) & (line.strip()[13:15] == 'CA') & (line.strip().startswith('ATOM')):
                cord2 = line.strip()[32:55].split()
#                     print(cord2)
            elif (line.strip()[22:27].strip() == str(residue_ID3).strip()) & (line.strip()[13:15] == 'CA') & (line.strip().startswith('ATOM')):
                cord3 = line.strip()[32:55].split()
#                     print(cord3)
            elif (line.strip()[22:27].strip() == str(residue_ID4).strip()) & (line.strip()[13:15] == 'CA') & (line.strip().startswith('ATOM')):
                cord4 = line.strip()[32:55].split()
#                     print(cord4)

        cord = list(map(float, cord1)) + list(map(float, cord2)) +list(map(float, cord3)) + list(map(float, cord4))
        arry_cord = np.array(cord).reshape(4,3)
        cord_aver = np.average(arry_cord,axis=0)
        box_center = list(cord_aver)
####################################################################################
def box_center_ligand(reference,folder):
    #if user uses ligand to specify the location of pocket
    cmd.load(folder+reference)
    box_center = cmd.centerofmass(reference.split('.')[0])
    cmd.delete('all')
    return box_center

#####################################################################################            
def big_box(rec_file, lig_file, residue_ID1, residue_ID2, residue_ID3, residue_ID4):
    global box_center
    
    with open(rec_file) as receptor: 
        # content = receptor.readlines()
        # for line in content:
        #     if (line.strip()[22:27].strip() == str(centerRes_ID).strip())  & (line.strip()[13:15] == 'CA'):
        #         pro_center = line.strip()[32:55].split()
                # box_center = list(map(float, pro_center))
        cmd.load(receptor)
        cmd.select('pocket', 'resid ' + residue_ID1.strip() + '+ ' + residue_ID2.strip() + '+ ' + residue_ID3.strip() + '+ ' + residue_ID4.strip())
        box_center = cmd.centerofmass(pocket)
        cmd.delete('all')
                    
#####################################################################################
def write_config(rec_file, lig_file,folder,box_center,configure_file, Rg='1.7'):
    #This function is for creating subdirectory and corresponding config.txt file for docking
    #global name_dir
    #get the path of current working directory
    #cwd = os.getcwd() + '/'
#create a sub directory for each receptor + ligand
    #name the sub directory as rec_lig; for example: 4agd_b49
    #name_dir = str(rec_file).split('.')[0] + '_' + str(lig_file).split('.')[0]
    # #make the sub directory
    # try:
    #     os.mkdir(name_dir)
    # #name the input configure file as rec_lig_config.txt; for example: 4agd_b49_config.txt
    # except FileExistsError:
    #     print('directory exists.')        
    # name_config = name_dir + '_config.txt'
    # #move the rec_file and lig_file to the sub directory
    # source = cwd
    # destination = cwd + name_dir + '/'
    # shutil.copy(source + rec_file, destination + rec_file)
    # shutil.copy(source + lig_file, destination + lig_file)
    #write the configure file in the sub directory
    with open(configure_file, 'w') as config:
        #the receptor file declare
        config.write('receptor = ' + folder+ str(rec_file).split('.')[0] +'.pdbqt' + '\n')
        #the ligand file declare
        config.write('ligand = ' + folder+ str(lig_file).split('.')[0] +'.pdbqt' + '\n')
        config.write('\n')
        #the docking box center declare
        config.write('center_x = ' + str(box_center[0]) +'\n')
        config.write('center_y = ' + str(box_center[1]) +'\n')
        config.write('center_z = ' + str(box_center[2]) +'\n')
        config.write('\n')
        #the docking box size declare: 
        #if no active site information is given, draw a big box
        # if (not (residue_ID1 or residue_ID2 or residue_ID3 or residue_ID4)) and (not reference):
        #     config.write('size_x = ' + '40' + '\n')
        #     config.write('size_y = ' + '40' + '\n')
        #     config.write('size_z = ' + '40' + '\n')
        #     config.write('\n')
        #     config.write('exhaustiveness = 10\n')
        # # or draw a smaller box in active site
        # else:
            #calculate the optimal box of the ligand
        box_cmmd = 'perl enzyme_evolver/workflow_src/binding_mode/docking/eBoxSize.pl ' + folder + str(lig_file).split('.')[0] +'.pdbqt'
        print(box_cmmd)
        optimal_size = sp.run(box_cmmd, shell = True, stdout = sp.PIPE).stdout.strip().decode("utf-8")
        print(optimal_size)
        print(type(optimal_size))
        print(Rg)
        print(type(Rg))
        config.write('size_x = ' + str(float(Rg) * float(optimal_size.strip())) + '\n')
        config.write('size_y = ' + str(float(Rg) * float(optimal_size.strip())) + '\n')
        config.write('size_z = ' + str(float(Rg) * float(optimal_size.strip())) + '\n')
        config.write('\n')
        config.write('exhaustiveness = 10\n')
        # config.write('num_modes = 20')
            
#####################################################################################            
def run_vina(configure_file,logfile,folder):
    command = 'vina --config ' + folder+ configure_file + ' --out ' + folder+ configure_file.split('.')[0]+ '_out.pdbqt ' +' --log ' + folder+ logfile
    try:
        sp.run(command,shell=True)
    except:
        pass
def ligand_preprocess(lig_file, folder, add_h=False,gen_3D=True):
    oname = str(lig_file).split('.')[0]+'.mol2'
    if (add_h) and (not gen_3D):
        #babel add hydrogen
        command_babel = f'babel {folder}{lig_file} {folder}{oname} -h -p 7'
        print(command_babel)
        sp.run(command_babel, shell=True)
    if (not add_h) and ( gen_3D):
        #babel not add hydrogens
        command_babel = f'babel {folder}{lig_file} {folder}{oname} --gen3D -p 7'
        print(command_babel)
        sp.run(command_babel,shell=True)
    if (add_h) and ( gen_3D):
        command_babel = f'babel {folder}{lig_file} {folder}{oname} --gen3D -h -p 7'
        print(command_babel)
        sp.run(command_babel, shell=True)
    if (not add_h) and ( not gen_3D):
        command_babel = f'babel {folder}{lig_file} {folder}{oname}'
        print(command_babel)
        sp.run(command_babel, shell=True)

def run_docking(rec_file, lig_file, folder, ref_ligand='', add_h=False,gen_3D=True, Rg='1.7'):

    #babel processing
    ligand_preprocess(lig_file=lig_file, folder=folder,add_h=add_h,gen_3D=gen_3D)
    #prepare ligands and receptors#
    prepare(rec_file, lig_file, folder)
    if ref_ligand == '':
        #find binding sites by epos
        find_site.site_calculate('aln_'+rec_file.strip().split('.')[0]+'.pdb',folder)
        #find the center of each site
        centers = find_site.centers_patch(folder)
    else:
        center = box_center_ligand(ref_ligand,folder)
        centers = []
        centers.append(center)

    print(centers)
    #dock ligand to each site
    for center in centers:
        #write configure file, output config1.txt coresponding to the the first binding site
        write_config(rec_file,lig_file,folder,center,folder+'config'+ str(centers.index(center)+1) + '.txt',Rg)
        #run vina docking
        run_vina('config'+ str(centers.index(center)+1) + '.txt', 'log'+ str(centers.index(center)+1) + '.txt',folder)

def run_docking_no_ali_prep(rec_file, lig_file, folder,ref_ligand='',add_h=False,gen_3D=True,Rg='1.7'):
    #prepare ligands
    ligand_preprocess(lig_file=lig_file, folder=folder,add_h=add_h,gen_3D=gen_3D)
    comand_lig_prepare = 'python2.5 enzyme_evolver/workflow_src/binding_mode/docking/prepare_ligand4.py -l ' + folder + \
                         str(lig_file).split('.')[0] + '.mol2' + ' -o ' \
                         + folder + str(lig_file).split('.')[0] + '.pdbqt'
    os.system(comand_lig_prepare)
    #no preparing proteins
    if ref_ligand == '':
        #find binding sites by epos
        find_site.site_calculate(rec_file.strip().split('.')[0]+'.pdb',folder)
        #find the center of each site
        centers = find_site.centers_patch(folder)
    else:
        center = box_center_ligand(ref_ligand,folder)
        centers = []
        centers.append(center)

    #print(centers)
    #dock ligand to each site
    for center in centers:
        #write configure file, output config1.txt coresponding to the the first binding site
        write_config(rec_file,lig_file,folder,center,folder+'config'+ str(centers.index(center)+1) + '.txt',Rg)
        #run vina docking
        run_vina('config'+ str(centers.index(center)+1) + '.txt', 'log'+ str(centers.index(center)+1) + '.txt',folder)

def rec_score(docking_folder):
    command = 'sh enzyme_evolver/workflow_src/binding_mode/docking/extract_score.sh ' + docking_folder
    try:
        sp.run(command,shell=True)
    except:
        pass
def summary_score(ligname,docking_folder,folder):
    command = 'head -1 ' + docking_folder+'scores.txt >>' + folder+ligname+'_summary_score.txt' + ';' \
                                                                                         'sort -n -k2 ' + folder+ligname+'_summary_score.txt > ' + folder +ligname+ '_ranking.txt'
    try:
        sp.run(command,shell=True)
    except:
        pass
def best_mode(rec_name,lig_name, docking_folder,folder):
    try:
        sp.run('sh enzyme_evolver/workflow_src/binding_mode/docking/extract_mode.sh ' + rec_name + ' ' + lig_name +' ' + docking_folder + ' '+ folder,shell=True)
        mode_name = f"{rec_name}_{lig_name}_best_mode.pdbqt"
        obj_name = f"{rec_name}_{lig_name}_best_mode"
        cmd.load(f"{folder}{mode_name}")
        cmd.split_states(f"{obj_name}")
        cmd.delete(f"{obj_name}")
        cmd.multisave(f"{folder}{obj_name}.pdb")
        cmd.delete('all')
        #sp.run(f"rm {folder}{mode_name}",shell=True)
        sp.run(f"sed -i '/HEADER/d' {folder}{obj_name}.pdb", shell=True)
    except:
        pass





