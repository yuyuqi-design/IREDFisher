from enzyme_evolver.workflow_main.evolutionary_analysis_main import EvolAnalysis
from enzyme_evolver.workflow_main.HomoModelling_main import HomoModelling
from enzyme_evolver.workflow_main.dockin_main import Docking
import subprocess as sp
import os
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
import shutil
from pathlib import Path

class Panel():
    def __init__(self,folder_id,protein_name, protein_seq,complex_struct='complex.pdb', cof_name='COF', lig_name='LIG',auto_single=True, homodimer=False,start_residue='', end_residue='', panel_sequence=False, test_mode=False, print_log=True):
        self.masterpath = Path.cwd()
        self.path = f"{Path(__file__).parents[1]}"
        self.folder_id = folder_id
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.protein_name = protein_name
        self.protein_seq = protein_seq
        self.protein_fasta = 'protein.fasta'
        self.ligands_file = 'ligands.txt'
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.complex_struct = complex_struct
        self.cof_name = cof_name.strip()
        self.lig_name = lig_name.strip()
        self.auto_single = auto_single
        self.homodimer = homodimer
        self.start_residue = start_residue
        self.end_residue = end_residue
        self.panel_sequence = panel_sequence
        self.db_job = Job.objects(folder_id=folder_id)[0]
        self.db_job.update_notes(f"panel_optimization = {panel_sequence}, modelling default = {auto_single}, homodimer modelling = {homodimer}")
        self.test_mode = test_mode
        self.print_log = print_log

        self._log(f"Folder ID = {self.folder_id}")
        self._log(f"Folder path = {self.folder_path}")

    def save_fasta(self):
        path = f"{self.folder_path}/{self.protein_fasta}"
        with open(path, 'w') as file:
            file.write('>')
            file.write(self.protein_name)
            file.write('\n')
            file.write(self.protein_seq)
            file.write('\n')
        sp.run(f'dos2unix {path}', shell=True)
        self.db_job.add_file(self.protein_fasta)

    def lig_ref_extract(self):
        #prepare the complex.pdb
        sp.run(f"sed -i 's/HETATM/ATOM  /g' {self.folder_path}/{self.complex_struct}", shell=True)
        #extract cof and lig
        with open(f"{self.folder_path}/{self.complex_struct}") as ref_complex:
            ref_complex_content = ref_complex.readlines()
            lig_lines = []
            cof_lines =[]
            for line in ref_complex_content:
                if line.startswith('ATOM'):
                    #print(self.lig_name)
                    if line[17:20] == self.lig_name:
                        #print(line)
                        lig_lines.append(line)
                    if line[17:20] == self.cof_name:
                        #print(line)
                        cof_lines.append(line)

            with open(f"{self.folder_path}/lig_ref.pdb", 'w') as lig:
                lig.writelines(lig_lines)
            with open(f"{self.folder_path}/cof_ref.pdb", 'w') as cof:
                cof.writelines(cof_lines)
        rm_lig_cof = f'sed -i "/{self.lig_name}/d;/{self.cof_name}/d" {self.folder_path}/{self.complex_struct}'
        print(rm_lig_cof)
        sp.run(rm_lig_cof, shell=True)

    def put_cof(self,rec_file,cof):
        with open(cof) as cof:
            coff_str = cof.read()

        f = open(rec_file, 'a')
        f.write(coff_str)
        f.close()



    def make_zip(self):
        self.db_job.update_status('Making zip file', print_log=self.print_log)
        os.chdir(f'{self.folder_path}')
        zipfile = self.folder_id +'.zip'
        sp.run(f'zip -r {zipfile} *', shell=True)
        os.chdir(self.masterpath)
        self.db_job.add_file(f'{self.folder_id}.zip')

    def job_finished(self):
        self.db_job.update_status('Job finished', print_log=self.print_log)

    def _log(self, msg):
        if self.print_log is True:
            print(msg)

    def run(self):
        if self.protein_seq:
            ## step 1
            self.db_job.update_status(f'finding homologous sequneces', print_log=self.print_log)
            evol_analysis = EvolAnalysis(self.folder_id, test_mode=self.test_mode, print_log=self.print_log, trim=True)
            self.save_fasta()
            evol_analysis.run()
            sp.run(f"cp {self.folder_path}/protein.fasta {self.folder_path}/allseq.fasta", shell=True)
            sp.run(f"cat {self.folder_path}/homologous_sequences75.fas >> {self.folder_path}/allseq.fasta", shell=True)
        else:
            sp.run(f"cp {self.folder_path}/{self.panel_sequence} {self.folder_path}/allseq.fasta", shell=True)
        ## step 2
        self.db_job.update_status(f'building 3D models for representative sequences', print_log=self.print_log)
        hm = HomoModelling(self.folder_id, auto_single=self.auto_single, auto_multiple=False, template_single=False, homodimer=self.homodimer,
                           start_residue=self.start_residue, end_residue=self.end_residue)
        fasta_string = open(f'{self.folder_path}/allseq.fasta', 'r').read()
        hm.create_sequence_files(fasta_string)
        hm.run()
        # # list the models.pdb > receptors.txt
        # rename the models' names


        ##step 3
        self.db_job.update_status(f'dock the substrate into active site of each homologs', print_log=self.print_log)
        with open(f'{self.folder_path}/receptors.txt', 'w') as f:
            f.write(self.complex_struct + '\n')
        sp.run(f"sed 's/$/\.pdb/g;$d' {self.folder_path}/codes.txt >> {self.folder_path}/receptors.txt", shell=True)
        if self.cof_name:
            docking_lig = Docking(self.folder_id, ref_ligand='lig_ref.pdb')
            receptors, ligands = docking_lig.get_rec_lig()
            docking_lig.align_rec(receptors)
            sp.run(f"sed -i '/END/d' {self.folder_path}/aln_*", shell=True)
            ###put cofactor
            for rec in receptors:
                aln_rec = f"{self.folder_path}/aln_{rec}"
                cof_file = f"{self.folder_path}/cof_ref.pdb"
                # print(aln_rec)
                # print(cof_file)

                self.put_cof(aln_rec, cof_file)
            #            print(receptors)
            #           print(ligands)
            docking_lig.docking_main(receptors, ligands)
        else:
            docking_lig = Docking(self.folder_id, ref_ligand='lig_ref.pdb')
            docking_lig.run()


        # Summariz
        self.db_job.update_status(f'Ranking sequences', print_log=self.print_log)
        self.make_zip()
        self.job_finished()







if __name__ == '__main__':
    from enzyme_evolver.mongo.default_connection import make_default_connection
    make_default_connection()

    from enzyme_evolver.app.workflow.functions import job_functions
    #folder_id = job_functions.create_new_job('test_protein', 'screen_panel', no_user=True)
    folder_id = '5a90ef28-1a62-4da3-922d-56833f677d4f'
    path_data = '/home/g02808yy/data/webserver/EnzymeEvolver3/EnzymeEvolver/enzyme_evolver/database'
    ###ligands.txt and ligand structure files
    sp.run(f"cp {path_data}/ligands.txt {path_data}/{folder_id}/", shell=True)
    sp.run(f"cp {path_data}/PCH.pdb {path_data}/{folder_id}/", shell=True)
    #evol_analysis = EvolAnalysis(folder_id, test_mode=True, print_log=True)
    #input sequence
    panel = Panel(folder_id,complex_struct='complex.pdb', cof_name='COF', lig_name='LIG')
    panel.save_fasta('protein', 'MSDQPRPTVIITGASSGVGLYATKALANRGWHVIMACRNLEKAEQAAKNLQIPPEAYTILHLDLSSLASVRGFVESFRA'
                                'LNRPLRALVCNAAVYYPLLKEPIYSVDGYEITVATNHLGHFLLINLLLEDLKNSPESDKRLVILGTVTANRKELGGKIPIPAPPDL'
                                'GNLEGFEKGFKKPIAMINGKPFKSGKAYKDSKLCNMLTARELHRRFHESTGIVFNSLYPGCVADTPLFRHHFPLFQKLFPLFQKKIT'
                                'GGYVSQELAGERVAMVVADPEFRQSGVHWSWGNRQKEGRKAFVQELSAEASDEQKARRLWELSEKLVGLA')
    panel.lig_ref_extract()
    #get string from form: cof_name,lig_name
    panel.run(targetLig='PCH.pdb')