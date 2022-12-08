import subprocess as sp
import os
import fnmatch
import shutil
from enzyme_evolver.workflow_src.evolutionary_analysis import get_homo
from enzyme_evolver.workflow_src.evolutionary_analysis.sequence_process import sequence_process
from enzyme_evolver.workflow_src.evolutionary_analysis.seqlogo import seqlogo
from enzyme_evolver.workflow_src.evolutionary_analysis import batch_replace
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
from pathlib import Path

class EvolAnalysis():

    def __init__(self, folder_id, test_mode=False, print_log=True, trim=False, ancestral=False):
        self.masterpath = Path.cwd()
        self.path = f"{Path(__file__).parents[1]}"
        self.folder_id = folder_id
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.trim = trim
        self.bootstrap_ancestral = ancestral
        self.db_job = Job.objects(folder_id=folder_id)[0]
        self.test_mode = test_mode
        self.print_log = print_log
        self.protein_fasta = "protein.fasta"

        #self.db_job.update_notes(f"Trim mode = {trim}, Bootstrap ancestral mode = {ancestral}")

        self._log(f"Folder ID = {self.folder_id}")
        self._log(f"Folder path = {self.folder_path}")

    def save_fasta(self, protein_name, protein_seq):
        path = f"{self.folder_path}/{self.protein_fasta}"
        with open(path, 'w') as file:
            file.write('>')
            file.write(protein_name)
            file.write('\n')
            file.write(protein_seq)
            file.write('\n')
        sp.run(f'dos2unix {path}', shell=True)
        self.db_job.add_file(self.protein_fasta)

    def run(self):
        self.blast()
        df_blast = self.parse_blast()
        df_sum = self.get_hit_summary(df_blast)
        self.cluster_seqs()
        df_sele = self.dump_redundant_sequences(df_sum)
        self.get_taxonomies(df_sele)
        self.muscle_alignment()
        # self.make_tree()
        self.make_zip()
        self.job_finished()

    def blast(self):
        self.db_job.update_status('Running blast', print_log=self.print_log)
        if self.test_mode is False:
            path_to_fasta = f"{self.folder_path}/protein.fasta"
            get_homo.run_blast_remote(path_to_fasta, refseqdb='refseq_protein', ofolder=f"{self.folder_path}/", max_seq='1000')
        else:
            database_functions.test_blastxml(self.folder_id)

        self.db_job.add_file('blast.xml')

    def parse_blast(self):
        self.db_job.update_status('Parsing blast', print_log=self.print_log)
        df_blast = get_homo.parse_blast(file=f"{self.folder_path}/blast.xml", save_csv=False)
        return df_blast

    def get_hit_summary(self, df_blast):
        self.db_job.update_status('Getting hit summary', print_log=self.print_log)
        df_sum = get_homo.hit_summary(df_blast, save_csv=False, ofolder=f"{self.folder_path}/")
        return df_sum

    def cluster_seqs(self):
        self.db_job.update_status('Clustering sequences by 90% similarity', print_log=self.print_log)
        get_homo.seqcluster(seq='homologous_sequences.fas', ofolder=f"{self.folder_path}/")
        self.db_job.add_file('cluster_organism.txt')

    def dump_redundant_sequences(self, df_sum):
        self.db_job.update_status('Dump redundant sequences', print_log=self.print_log)
        df_sele = get_homo.sequence_dump(df=df_sum, seq_dump=f"{self.folder_path}/cluster_organism.txt")
        return df_sele

    def get_taxonomies(self, df_sele):
        self.db_job.update_status('Get taxonomy of each organism', print_log=self.print_log)

        # get taxonomy of each organism
        taxoNCBI = f"{self.path}/database/taxoNCBI.csv"

        df_seq_tax = get_homo.taxo(df=df_sele, save_csv=False, ofsummary=f'{self.folder_path}/hit_summary.csv',
                                   ofseq=f'{self.folder_path}/homologous_sequences75.fas',
                                   taxoNCBI=taxoNCBI)

        # Hack to delete first line of homologous_sequences75.fas.
        with open(f'{self.folder_path}/homologous_sequences75.fas', 'r') as fin:
            data = fin.read().splitlines(True)
        if data[0] == 'Sequence':
            with open(f'{self.folder_path}/homologous_sequences75.fas', 'w') as fout:
                fout.writelines(data[1:])

        spec_char_script = f"{self.path}/workflow_src/evolutionary_analysis/rm_specialChara.sh"
        sp.run(f'sh {spec_char_script} {self.folder_path}/hit_summary.csv', shell=True)
        sp.run(f'sh {spec_char_script} {self.folder_path}/homologous_sequences75.fas', shell=True)

        self.db_job.add_file('hit_summary.csv')
        self.db_job.add_file('homologous_sequences75.fas')

        # creates a file for rename the phylogenetic tree tips
        batch_replace.replace_extract(df_seq_tax, f"{self.folder_path}/")
        self.db_job.add_file('replace.txt')

    def muscle_alignment(self):
        self.db_job.update_status('Multiple sequence alignment by muscle', print_log=self.print_log)

        sequence_process(ofolder=f"{self.folder_path}/", file='homologous_sequences75.fas', trim=self.trim,
                         ofseq_ali='homologous_sequences75_muscle.fas', ofhtml='homologous_sequences75_muscle.html')
        seqlogo(file='homologous_sequences75_muscle.fas', ofolder=f"{self.folder_path}/")

        self.db_job.add_file('sequence_logo.png')
        self.db_job.add_file('homologous_sequences75_muscle.fas')
        self.db_job.add_file('homologous_sequences75_muscle.html')

    def make_tree(self):

        self.db_job.update_status('Making phylogenetic tree', print_log=self.print_log)

        if self.trim == True:
            muscle_file = "homologous_sequences75_muscle.fas-gb"
        else:
            muscle_file = "homologous_sequences75_muscle.fas"

        if self.bootstrap_ancestral == False:
            tree_cmd = f"{self.path}/workflow_src/evolutionary_analysis/phylo.sh"
            tree_mode = "*support"
        else:
            tree_cmd = f"{self.path}/workflow_src/evolutionary_analysis/phylo_bs_ac.sh"
            tree_mode = "*bestTree"

        sp.run(f"sh {tree_cmd} {self.folder_path}/{muscle_file} {self.folder_path}/", shell=True)

        for treefile in os.listdir(self.folder_path):
            if fnmatch.fnmatch(treefile, f'{tree_mode}'):
                batch_replace.batch_replace(f"{self.folder_path}/{treefile}", f"{self.folder_path}/replace.txt")


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
    folder_id = job_functions.create_new_job('test_protein', 'Evolutionary Analysis', no_user=True)

    # test_mode requires a blast.xml file in /database
    evol_analysis = EvolAnalysis(folder_id, test_mode=True, print_log=True)
    evol_analysis.save_fasta('test_protein', 'MG')
    evol_analysis.run()
