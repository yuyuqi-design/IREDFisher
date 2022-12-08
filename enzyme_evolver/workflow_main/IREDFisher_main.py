from enzyme_evolver.workflow_main.evolutionary_analysis_main import EvolAnalysis
from enzyme_evolver.workflow_main.HomoModelling_main import HomoModelling
from enzyme_evolver.workflow_main.dockin_main import Docking
from enzyme_evolver.workflow_src.binding_mode.docking import autodocking
from enzyme_evolver.workflow_src.IREDFisher_scoring import IREDFisher_scoring
import subprocess as sp
import os
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
import shutil
from pathlib import Path

class IREDFisher():
    def __init__(self,folder_id,database = 'Screened', N_index = '3', panel_sequence=False, complex_struct='complex.pdb', start_residue='9', end_residue='294',  test_mode=False, print_log=True):
        self.masterpath = Path.cwd()
        self.path = f"{Path(__file__).parents[1]}"
        self.folder_id = folder_id
        self.database = database
        self.N_index = N_index
        self.panel_sequence = panel_sequence
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.ligands_file = 'ligands.txt'
        self.complex_struct = complex_struct
        self.start_residue = start_residue
        self.end_residue = end_residue
        self.db_job = Job.objects(folder_id=folder_id)[0]
        self.db_job.update_notes(f"Fishing IREDs from {database} database, panel_optimization = {panel_sequence}")
        self.test_mode = test_mode
        self.print_log = print_log
        self._log(f"Folder ID = {self.folder_id}")
        self._log(f"Folder path = {self.folder_path}")


    def run_modelling(self):
        self.db_job.update_status(f'building 3D models for panel sequences', print_log=self.print_log)

        hm = HomoModelling(self.folder_id, auto_single=False, auto_multiple=False, template_single=False,
                           homodimer=True,
                           start_residue=self.start_residue, end_residue=self.end_residue)
        fasta_string = open(f'{self.folder_path}/allseq.fasta', 'r').read()
        hm.create_sequence_files(fasta_string)
        hm.run(IREDFisher=True)
    def prepare_rec_file(self):
        self.db_job.update_status(f'preparing the receptors.txt file', print_log=self.print_log)
        sp.run(f"cp {self.path}/database/complex.pdb {self.path}/database/complex.fasta {self.path}/database/lig_ref.pdb {self.path}/database/cof_ref.pdb {self.folder_path}", shell=True)
        with open(f'{self.folder_path}/receptors.txt', 'w') as f:
            f.write(self.complex_struct + '\n')
        sp.run(f"sed 's/$/\.pdb/g;/^$/d' {self.folder_path}/codes.txt >> {self.folder_path}/receptors.txt", shell=True)


    def run_docking(self):

        self.db_job.update_status(f'dock the substrate into active site of each homologs', print_log=self.print_log)
        docking_lig = Docking(self.folder_id, ref_ligand='lig_ref.pdb', cof = 'cof_ref.pdb')
        if self.panel_sequence:
            docking_lig.run()
        ###if use database, will use prepared models, no not need to prepare again.
        else:
            docking_lig.run_no_aln_prep()

    def run_scoring(self):
        self.db_job.update_status(f'Scoring sequences', print_log=self.print_log)
        Rescore = IREDFisher_scoring(self.folder_id)
        Rescore.run_rescore()
        self.db_job.update_status(f'Scoring finished', print_log=self.print_log)

        ########

    def job_finished(self):
        self.db_job.update_status('Job finished', print_log=self.print_log)

    def _log(self, msg):
        if self.print_log is True:
            print(msg)
    def make_zip(self):
        self.db_job.update_status('Making zip file', print_log=self.print_log)
        os.chdir(f'{self.folder_path}')
        zipfile = self.folder_id + '.zip'
        sp.run(f'zip -r {zipfile} *', shell=True)
        os.chdir(self.masterpath)
        self.db_job.add_file(f'{self.folder_id}.zip')

    def run_IREDFisher(self):

        if self.panel_sequence:
            sp.run(f"cp {self.path}/database/singleT.pdb {self.folder_path}", shell=True)
            sp.run(f"cp {self.path}/database/dimerT.pdb {self.folder_path}", shell=True)
            self.run_modelling()
            self.prepare_rec_file()
        else:
            sp.run(f"cp {self.path}/database/{self.database}/* {self.folder_path}", shell=True)
            #receptors named without aln_
            ###use aligned prepared models(pdbqt) with NAP bound, receptors (1.pdbqt 2.pdbqt, 3.pdbqt)

        sp.run(
            f"cp {self.path}/database/complex.pdb {self.path}/database/lig_ref.pdb {self.path}/database/cof_ref.pdb {self.folder_path}",
            shell=True)
        self.run_docking()
        self.run_scoring()
        #remove folders
        sp.run(f"rm -rf {self.folder_path}/*/", shell=True)
        self.make_zip()
        self.job_finished()


