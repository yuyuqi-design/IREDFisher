import subprocess as sp
import sys
import parmed as pmd


def lig_process(file,folder,lig_charge='0'): #file pdb file
    with open(folder+file) as lig:
        file = file.strip()
        with open(folder+'lig_H.pdb', 'w') as out:
            content = lig.readlines()
            ##replace HETATM with ATOM
            if 'HETATM' in content[3]:
                atoms = [lines for lines in content if lines.startswith('HETATM')]
                for atom in atoms:
                    atom = 'ATOM  ' + atom[6:]
                    out.write(atom)
            else:
                out.writelines(content)

        #add_H = 'babel -h ' + folder + file +'_out.pdb ' + folder+'lig_H.pdb -p 7.4'
        #sp.run(add_H, shell=True)

        ff = 'antechamber -i '+ folder + 'lig_H.pdb' +' -fi pdb -o '+ folder + 'lig.mol2' +' -fo mol2 -c bcc -s 2 -nc ' +lig_charge
        sp.run(ff, shell=True)
        parmchk = 'parmchk2 -i '+ folder + 'lig.mol2' +' -f mol2 -o '+ folder + 'lig.frcmod'
        sp.run(parmchk, shell=True)
        del_conect = "sed -i '/CONECT/d' " + folder + 'lig_H.pdb'
        sp.run(del_conect,shell=True)
        #sp.run('rm lig_out.pdb',shell=True)
        sp.run('rm ANTECHAMBER_*',shell=True)

def assemble(rec_file,lig_file='lig_H.pdb',folder=''):
    with open(folder+rec_file) as pro:
        with open(folder+lig_file) as lig:
            with open(folder+'complex.pdb', 'w') as output:
                pro_content = pro.readlines()
                lig_content = lig.readlines()
                complex_content = pro_content + lig_content
                output.writelines(complex_content)
    sp.run(f'sed -i "/END/d" {folder}complex.pdb',shell=True)


def ff_tleap(file='complex.pdb', ligname='MOL',folder = ''): #give a complete structure as file input (all protein chains,ligand)
    file = file.strip()
    f_content = open('enzyme_evolver/workflow_src/binding_mode/md/tleap.in').read()
    f_content = f_content.replace('ligname' , ligname.strip())
    f_content = f_content.replace('folder',folder.strip())
    with open(folder+'tleap2.in','w') as replace:
        replace.writelines(f_content)
    if file != 'complex.pdb':
        sp.run('cp ' + folder+file.strip() + ' '+ folder+'complex.pdb', shell = True)
    tleap = 'tleap -f ' + folder+ 'tleap2.in'
    sp.run(tleap, shell=True)


def ff_tleap_apo(file='pro.pdb', folder = ''):
    file = file.strip()
    f_content = open('enzyme_evolver/workflow_src/binding_mode/md/tleap_apo.in').read()
    f_content = f_content.replace('folder',folder.strip())
    with open(folder+'tleap2.in','w') as replace:
        replace.writelines(f_content)
    if file != 'pro.pdb':
        sp.run('cp ' + folder+file.strip() + ' '+ folder+'pro.pdb', shell = True)
    tleap = 'tleap -f ' + folder+ 'tleap2.in'
    sp.run(tleap, shell=True)



def amber2gmx(kind='complex',folder=''):
    amber = pmd.load_file(folder+kind+'.prmtop', folder+kind+'.inpcrd')
    # Save a GROMACS topology and GRO file
    try:
        amber.save(folder+kind+'.top')
        amber.save(folder+kind+'.gro')
    except IOError:
        pass


def gmx_restraint(input='complex.gro', folder=''):
    # command

    restr_Protein_comm = 'echo 1|gmx genrestr -f ' + str(folder + input) + ' -o ' + folder + 'posre_PROTEIN.itp'
    restr_Heavy_comm = 'echo 2|gmx genrestr -f ' + str(folder + input) + ' -o ' + folder + 'posre_Heavy.itp'
    restr_MC_comm = 'echo 5|gmx genrestr -f ' + str(folder + input) + ' -o ' + folder + 'posre_MC.itp'
    restr_CA_comm = 'echo 3|gmx genrestr -f ' + str(folder + input) + ' -o ' + folder + 'posre_CA.itp'

    # call the commmand genrestr to generate restraint file

    print('starting to generate protein restraint...\n')
    sp.run(restr_Protein_comm, shell=True)
    print('starting to generate Heavy atoms restraint...\n')
    sp.run(restr_Heavy_comm, shell=True)
    print('starting to generate MC atoms restraint...\n')
    sp.run(restr_MC_comm, shell=True)
    print('starting to generate CA atoms restraint...\n')
    sp.run(restr_CA_comm, shell=True)
    sp.run('rm \#*', shell=True)

def lig_restraint(input='lig_H.pdb', folder=''):
    restr_lig_comm = 'echo 1|gmx genrestr -f ' + str(folder + input) + ' -o ' + folder + 'posre_LIG.itp'
    sp.run(restr_lig_comm, shell=True)

def top_edit(input='complex.top',folder='', apo=True):
    with open(f'{folder}/{input}') as fin:
        with open(f'{folder}/topol_restraint.top', 'w') as fout:
            top0 = fin.readlines()
            print(top0)
            ##top molecules order : system1 (protein), ligand, ions, water
            ##get the line number of 'moleculetyle'
            molecule_line = []
            for i, x in enumerate(top0):
                if x == '[ moleculetype ]\n':
                    molecule_line.append(i)
            #print(molecule_line)
            #insert the text for restrain protein
            protein_restraint = f'; Include Position restraint file\n#ifdef POSRES_PROTEIN\n' \
                                f'#include "posre_PROTEIN.itp"\n#endif\n#ifdef POSRES_Heavy\n' \
                                f'#include "posre_Heavy.itp"\n#endif\n#ifdef POSRES_MC\n' \
                                f'#include "posre_MC.itp"\n#endif\n#ifdef POSRES_CA\n#include "posre_CA.itp"\n' \
                                f'#endif\n\n'
            print(molecule_line[1])
            top0.insert(molecule_line[1], protein_restraint)
            ##if it is a ligand-bound complex, insert the ligand restraint
            if apo==False:
                ligand_restraint = f'; Include Position restraint file\n#ifdef POSRES_LIG\n' \
                                   f'#include "posre_LIG.itp"\n#endif\n\n'
                top0.insert(molecule_line[2], ligand_restraint)
            fout.writelines(top0)
if __name__ == '__main__':
    top_edit(input='complex.top',folder='/home/g02808yy/data/webserver/Github/IREDFisher/enzyme_evolver/database/7de1f30b-a284-4f30-9e9f-21be48a6f587/5mst_AMP/', apo=False)







        

