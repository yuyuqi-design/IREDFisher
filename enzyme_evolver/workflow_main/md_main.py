from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
from enzyme_evolver.workflow_src.binding_mode.docking.mkdir import mk_dir
from enzyme_evolver.workflow_src.binding_mode.md import prepare_md
import sys
import subprocess as sp
from pathlib import Path
import os

#cat file.txt
#6h44.pdb apo (only protein system)
#or
#6h44.pdb Trp.pdb (protein-ligand system)
#protein hydrogens needs to be deleted and ligand should be a hydrogen-added molecule
#(in the case of ASH, you need to specify the residue name to indicate the protonation state.)

class MDpreparation():
    def __init__(self, folder_id, job_detail='file.txt', test_mode=False, print_log=True):
        # masterpath is under /EnzymeEvolver
        self.masterpath = Path.cwd()
        # path is under /EnzymeEvolver/enzyme_evolver
        self.path = f"{Path(__file__).parents[1]}"
        self.folder_id = folder_id
        # folder_path is each jobs's folder
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.db_job = Job.objects(folder_id=folder_id)[0]
        self.test_mode = test_mode
        self.print_log = print_log
        self.job_detail = 'file.txt'

        self._log(f"Folder ID = {self.folder_id}")
        self._log(f"Folder path = {self.folder_path}")


    def run(self):
        folder = f'{self.folder_path}/'
        with open(folder + 'file.txt') as struc_file:
            structures = struc_file.readlines()
            for structure in structures:
                # cat file.txt
                # 6h44.pdb apo (only protein system)
                # or
                # 6h44.pdb Trp.pdb (protein-ligand system)
                rec_file = structure.split()[0].strip()
                lig_file = structure.split()[1].strip()
                lig_name = structure.split()[2].strip()
                lig_charge = structure.split()[3].strip()
                # print(rec_file + lig_file)
                # docking_main(receptors,ligands,folder,ref_ligand='lig_ref.pdb')
                self.md_main(rec_file, lig_file, ligname=lig_name, lig_charge=lig_charge)
        self.make_zip()
        self.job_finished()

    def md_main(self, rec_file, lig_file, ligname = 'MOL',lig_charge = '0', cof1 = '', cof1_charge = '', cof2 = '', cof2_charge = ''):
        folder = f'{self.folder_path}/'

        #print(rec_file)
        #print(lig_file)


    #############################################################################

        if lig_file == 'apo':#only protein system

            name_dir = mk_dir(rec_file,'apo.pdb',folder)
            print(name_dir)
            sub_folder = folder + name_dir + '/'
            sp.run('cp ' + folder + rec_file + ' '+ folder + name_dir, shell=True)
            sp.run(f'cp {self.path}/workflow_src/binding_mode/md/job0.sh {sub_folder}', shell=True)
            sp.run(f'cp {self.path}/workflow_src/binding_mode/md/*mdp {sub_folder}', shell=True)
            self.db_job.update_status('running tleap', print_log=self.print_log)
            prepare_md.ff_tleap_apo(rec_file, folder = sub_folder)
            self.db_job.update_status('parameters are being converted to gromacs format', print_log=self.print_log)
            prepare_md.amber2gmx(kind='pro',folder=sub_folder)
            self.db_job.update_status('putting restraints in parameters files', print_log=self.print_log)
            prepare_md.gmx_restraint(input='pro.gro',folder=sub_folder)
            self.db_job.update_status('editing the topology file', print_log=self.print_log)
            prepare_md.top_edit(input='pro.top',folder=sub_folder,apo=True)
            commd_itp = f"sed -i 's/-DPOSRES_LIG//g' {sub_folder}*.mdp"
            sp.run(commd_itp, shell=True)


        else:#protein-ligand system
            name_dir = mk_dir(rec_file,lig_file,folder)
            print(name_dir)
            sub_folder = folder + name_dir + '/'

            sp.run(
                'cp ' + folder + lig_file + ' ' + folder +  rec_file + ' ' + folder + name_dir, shell=True)
            sp.run(f'cp {self.path}/workflow_src/binding_mode/md/job0.sh {sub_folder}', shell=True)
            sp.run(f'cp {self.path}/workflow_src/binding_mode/md/*mdp {sub_folder}', shell=True)
            #process the ligand, add hydrogens, generating AM1-BCC charges (lig.mol2) and force field files (lig.frcmod) and a pdb lig_H.pdb
            self.db_job.update_status('calculating parameter for ligand using antechamber', print_log=self.print_log)
            prepare_md.lig_process(file=lig_file,folder=sub_folder,lig_charge=lig_charge)
            #output complex.pdb
            self.db_job.update_status('putting ligand and protein together', print_log=self.print_log)
            prepare_md.assemble(rec_file=rec_file,lig_file='lig_H.pdb',folder=sub_folder)
            self.db_job.update_status('running tleap for the complex structure', print_log=self.print_log)
            prepare_md.ff_tleap(file='complex.pdb',ligname=ligname,folder = sub_folder)
            prepare_md.amber2gmx(folder=sub_folder)
            self.db_job.update_status('putting restraints in parameters files', print_log=self.print_log)
            prepare_md.gmx_restraint(input='complex.gro',folder=sub_folder)
            prepare_md.lig_restraint(input='lig_H.pdb',folder=sub_folder)
            self.db_job.update_status('editing the topology file', print_log=self.print_log)
            prepare_md.top_edit(input='complex.top', folder=sub_folder, apo=False)

    def make_zip(self):
        self.db_job.update_status('Making zip file', print_log=self.print_log)
        os.chdir(f'{self.folder_path}')
        zipfile = self.folder_id + '.zip'
        sp.run(f'zip -r {zipfile} *', shell=True)
        os.chdir(self.masterpath)
        self.db_job.add_file(f'{self.folder_id}.zip')

    def job_finished(self):
        self.db_job.update_status('Job finished', print_log=self.print_log)

    def _log(self, msg):
        if self.print_log is True:
            print(msg)



if __name__ == '__main__':
    fileid = 'dimer_thal'
    folder = 'enzyme_evolver/database/' + fileid+'/'
    with open(folder+'file.txt') as struc_file:
            structures = struc_file.readlines()
            for structure in structures:
                rec_file = structure.split()[0].strip()
                lig_file = structure.split()[1].strip()
                # print(rec_file + lig_file)
                #docking_main(receptors,ligands,folder,ref_ligand='lig_ref.pdb')
                md_main(rec_file, lig_file, folder, ligname = 'MOL',lig_charge = '0')


