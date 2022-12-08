from enzyme_evolver.workflow_src.binding_mode.docking import autodocking
import subprocess as sp
import shutil
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
from pathlib import Path
import os
#input a list; receptors, ligands which contains their corresponding names
class Docking:
    def __init__(self, folder_id, ref_ligand='', cof = '', add_h=False,gen_3D=False, Rg=1.7, test_mode=False, print_log=True):
        #masterpath is under /EnzymeEvolver
        self.masterpath = Path.cwd()
        #path is under /EnzymeEvolver/enzyme_evolver/
        self.path = f"{Path(__file__).parents[1]}"
        self.folder_id = folder_id
        #folder_path is each jobs's folder
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.db_job = Job.objects(folder_id=folder_id)[0]
        self.test_mode = test_mode
        self.print_log = print_log
        self.ref_ligand = ref_ligand
        self.cof = cof
        self.add_h = add_h
        self.gen_3D=gen_3D
        self.Rg = Rg


        self._log(f"Folder ID = {self.folder_id}")
        self._log(f"Folder path = {self.folder_path}")

    def get_rec_lig(self):
        self.db_job.update_status('getting the names of all receptors and ligands', print_log=self.print_log)
        #the receptors.txt and ligands.txt are uploaded by users along with the corresponding files.
        folder = f"{self.folder_path}/"
        sp.run(f'dos2unix {self.folder_path}/receptors.txt', shell=True)
        sp.run(f'dos2unix {self.folder_path}/ligands.txt', shell=True)
        with open(folder+'receptors.txt') as recs_file:
            with open(folder+'ligands.txt') as ligs_file:
                receptors = recs_file.readlines()
                receptors = [x.strip() for x in receptors]
                ligands = ligs_file.readlines()
                ligands = [l.strip() for l in ligands]
                print(receptors)
                print(ligands)
                return receptors,ligands

    def align_rec(self,receptors):
        self.db_job.update_status('aligning all receptors to the first protein', print_log=self.print_log)
        folder = f"{self.folder_path}/"
        autodocking.align(receptors, folder)
        for pdb in receptors:
            self.db_job.add_file('aln_'+pdb.strip().split('.')[0]+'.pdb')


    def put_cof(self,rec_file,cof):
        with open(cof) as cof:
            coff_str = cof.read()

        f = open(rec_file, 'a')
        f.write(coff_str)
        f.close()




    def docking_main(self, receptors, ligands, no_prepare=False):
        self.db_job.update_status('cross docking going on', print_log=self.print_log)
        for rec_file in receptors:
            for lig_file in ligands:
                #print(rec_file)
                #print(lig_file)
                folder = f"{self.folder_path}/"
                name_dir = autodocking.mk_dir(rec_file,lig_file,folder)
                sub_folder = folder + name_dir + '/'
                print(name_dir)
                try:
                    sp.run(f"cp {folder}{self.ref_ligand} {folder}{name_dir}", shell=True)
                except:
                    pass
                ##Rg value
                if self.Rg == '':
                    self.Rg = 1.7
                if no_prepare:
                    sp.run('cp ' + folder + lig_file + ' ' + folder  + rec_file + ' ' + folder + name_dir,
                           shell=True)
                    autodocking.run_docking_no_ali_prep(rec_file, lig_file, sub_folder, self.ref_ligand,add_h=self.add_h,gen_3D=self.gen_3D, Rg=self.Rg)
                else:
                    sp.run('cp ' + folder + lig_file + ' ' + folder + 'aln_' + str(rec_file).split('.')[0]+'.pdb' + ' ' + folder + name_dir,
                           shell=True)
                    autodocking.run_docking(rec_file,lig_file, sub_folder,self.ref_ligand,add_h=self.add_h,gen_3D=self.gen_3D,Rg=self.Rg)

                #calculate the docking score from each predicted binding site
                autodocking.rec_score(sub_folder)
                rec_name = str(rec_file).split('.')[0]
                lig_name = str(lig_file).split('.')[0]
                #extract the best binding mode and output to the folder directory
                autodocking.best_mode(rec_name,lig_name,sub_folder,folder)
                self.db_job.add_file(rec_name+'_'+lig_name+'_best_mode.pdb')
                #give a ranking score file
                autodocking.summary_score(lig_name,sub_folder,folder)
                # clean the directory
                #sp.run(f"rm -r {sub_folder}; {self.folder_path}/aln_*.pdb", shell=True)
            for lig_file in ligands:
                lig_name = str(lig_file).split('.')[0]
                self.db_job.add_file(lig_name + '_summary_score.txt')
                self.db_job.add_file(lig_name + '_ranking.txt')


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

    def run(self):
        receptors, ligands = self.get_rec_lig()
        self.align_rec(receptors)
        sp.run(f"sed -i '/END/d' {self.folder_path}/aln_*", shell=True)
        if self.cof:
            for rec in receptors:
                aln_rec = f"{self.folder_path}/aln_{rec}"
                cof_file = f"{self.folder_path}/cof_ref.pdb"
                sp.run(f"sed -i '/CONECT/d' {cof_file}", shell=True)
                # print(aln_rec)
                # print(cof_file)
                self.put_cof(aln_rec, cof_file)
        self.docking_main(receptors, ligands)
        self.make_zip()
        self.job_finished()
    def run_no_aln_prep(self):
        receptors, ligands = self.get_rec_lig()
        self.docking_main(receptors, ligands, no_prepare=True)
        self.make_zip()
        self.job_finished()




if __name__ == '__main__':
    # fileid = 'POR/docking'
    # folder = 'enzyme_evolver/database/' + fileid+'/'
    # with open(folder+'receptors.txt') as recs_file:
    #     with open(folder+'ligands.txt') as ligs_file:
    #         receptors = recs_file.readlines()
    #         receptors = [x.strip() for x in receptors]
    #         print(receptors)
    #         ligands = ligs_file.readlines()
    #         ligands = [l.strip() for l in ligands]
    #         print(ligands)
    #         docking_main(receptors,ligands,folder,ref_ligand='lig_ref.pdb')
            #docking_main(receptors,ligands,folder,ref_ligand=False)
    from enzyme_evolver.mongo.default_connection import make_default_connection
    make_default_connection()
    from enzyme_evolver.app.workflow.functions import job_functions

    folder_id = job_functions.create_new_job('test_protein', 'Docking', no_user=True)
    # folder_id = 'test'
    path_data = '/home/g02808yy/data/webserver/EnzymeEvolver3/EnzymeEvolver/enzyme_evolver/database'
    sp.run(f"cp {path_data}/receptors.txt {path_data}/{folder_id}/",shell=True)
    sp.run(f"cp {path_data}/ligands.txt {path_data}/{folder_id}/", shell=True)
    sp.run(f"cp {path_data}/*pdb {path_data}/{folder_id}/", shell=True)
    docking = Docking(folder_id)
    #hm.create_sequence_files(fasta_string)
    docking.run()



