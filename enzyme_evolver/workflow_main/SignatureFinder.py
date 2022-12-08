from enzyme_evolver.workflow_main.phylogeny_main import Phylogeny
from enzyme_evolver.workflow_main.HomoModelling_main import HomoModelling
from enzyme_evolver.workflow_main.dockin_main import Docking
from enzyme_evolver.workflow_src.binding_mode.docking import autodocking
from enzyme_evolver.workflow_src.binding_mode.docking.complex_split import complex_split
from enzyme_evolver.workflow_src.homology_modelling.rmsdCA_pro import split_positive_negative
from enzyme_evolver.workflow_src.evolutionary_analysis.sequence_process import sequence_process
from enzyme_evolver.workflow_src.homology_modelling.merge import merge
from enzyme_evolver.workflow_src.evolutionary_analysis.seqlogo import seqlogo
import subprocess as sp
import os
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
import shutil
from pathlib import Path

class SigFinder():
    def __init__(self,folder_id, sequences='panel_sequence.fasta', complex_struct='complex.pdb', Rg=0.7, rmsd=2,test_mode=False, print_log=True):
        self.masterpath = Path.cwd()
        self.path = f"{Path(__file__).parents[1]}"
        self.folder_id = folder_id
        self.panel_sequence = sequences
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.ligands_file = 'ligands.txt'
        self.complex_struct = complex_struct
        self.db_job = Job.objects(folder_id=folder_id)[0]
        self.db_job.update_notes(f"Finding sequence signatures based on protein-ligand interaction")
        self.test_mode = test_mode
        self.Rg = float(Rg)
        self.rmsd = float(rmsd)
        self.print_log = print_log
        self._log(f"Folder ID = {self.folder_id}")
        self._log(f"Folder path = {self.folder_path}")

    def run_SignatureFinder(self):
        self.split_lig()
        ##sequence alignment for the uploaded sequence
        self.run_phylogeny()
        self.run_modelling()
        self.prepare_rec_file()
        ##sequence alignment
        positive, negative = self.rmsd_homologs2ref()
        df = merge(f"{self.folder_path}/models_rmsd2ref.csv",f"{self.folder_path}/Modelling_template_summary.txt")
        df.to_csv(f"{self.folder_path}/structual_annotation.csv", index=False)
        positive.iloc[:, 0].to_csv(f"{self.folder_path}/receptors.txt", index=False, header=False)
        ###if rmsd < value, add the seqeunce to putative seq and then sequence alignment
        #otherwise, add the sequence to negative guess and then sequence alignment
        self.muscle_alignment(file='positive')
        self.muscle_alignment(file='negative')
        # self.run_docking()
        # ##analyze the docking pose, set scaling factor as 0.5 times of gyration radius for Corrin ring, get rmsd for each pose.
        # self.rmsd_lig()
        # self.job_finished()
        self.make_zip()

    def split_lig(self):
        complex_split(complex=f'{self.folder_path}/{self.complex_struct}',pro=f'{self.folder_path}/pro.pdb',lig=f'{self.folder_path}/lig.pdb')

    def run_phylogeny(self):
        sp.run(f'cp {self.folder_path}/allseq.fasta {self.folder_path}/input.fasta', shell=True)
        phylogeny = Phylogeny(folder_id=self.folder_id, alignment='sequence',tree_number = '10')
        phylogeny.run()

    def run_modelling(self):
        self.db_job.update_status(f'building 3D models for panel sequences', print_log=self.print_log)

        hm = HomoModelling(self.folder_id, auto_single=True, auto_multiple=False, template_single=False,
                           homodimer=False)
        fasta_string = open(f'{self.folder_path}/allseq.fasta', 'r').read()
        hm.create_sequence_files(fasta_string)
        hm.run(IREDFisher=False)

    def rmsd_homologs2ref(self):
        with open(f"{self.folder_path}/receptors.txt") as f:
            recs = f.readlines()
            positive, negative = split_positive_negative(recs=recs,folder=self.folder_path,rmsd=self.rmsd)
            return  positive, negative

    def muscle_alignment(self,file):
        self.db_job.update_status('Multiple sequence alignment by muscle', print_log=self.print_log)

        sequence_process(ofolder=f"{self.folder_path}/", file=file+'.fasta', trim=False,
                         ofseq_ali=file+'_muscle.fasta', ofhtml=file+'_muscle.html')
        seqlogo(file=file+'_muscle.fasta', ofolder=f"{self.folder_path}/")

        self.db_job.add_file(file + '_muscle_sequence_logo.png')
        self.db_job.add_file(file+'_muscle.fasta')
        self.db_job.add_file(file+'_muscle.html')


    def prepare_rec_file(self):
        self.db_job.update_status(f'preparing the receptors.txt file', print_log=self.print_log)
        with open(f'{self.folder_path}/receptors.txt', 'w') as f:
            f.write('pro.pdb\n')
            for file in os.listdir(self.folder_path):
                if (file.split('.')[1] == 'pdb') and (file.split('.')[0] !='lig') and (file.split('.')[0] !='pro') and (file.split('.')[0] !='complex'):
                    f.write(file + '\n')
        #

        #
        # list_receptor_cmd = f'ls -1 *pdb > {self.folder_path}/receptors.txt'
        # print(list_receptor_cmd)
        # sp.run(list_receptor_cmd,shell=True)
        # # with open(f'{self.folder_path}/receptors.txt', 'w') as f:
        # #     f.write('pro.pdb\n')
        # receptor_file = f"sed -i '/lig/d;/^$/d' {self.folder_path}/receptors.txt"
        # print(receptor_file)
        # sp.run(receptor_file, shell=True)
        with open(f'{self.folder_path}/ligands.txt', 'w') as f2:
            f2.write('lig.pdb\n')

    def run_docking(self):

        self.db_job.update_status(f'dock the ligand into active site of each homologs', print_log=self.print_log)
        docking_lig = Docking(self.folder_id, ref_ligand='lig.pdb', Rg=self.Rg)
        docking_lig.run()

    def rmsd_lig(self):
        pass

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
    from enzyme_evolver.mongo.default_connection import make_default_connection

    make_default_connection()

    from enzyme_evolver.app.workflow.functions import job_functions

    folder_id = '5f25bd0b-593b-4406-9081-2501ae5996f9'
    signature = SigFinder(folder_id=folder_id, sequences='allseq.fasta', complex_struct='complex.pdb', Rg=1.0, rmsd=2,test_mode=False, print_log=True)
    signature.run_SignatureFinder()


