############################a upload portal for user to use their own trimmed sequence alignment file.
import os
import subprocess as sp
import sys
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
import uuid
import fnmatch
from enzyme_evolver.workflow_src.evolutionary_analysis import batch_replace
from enzyme_evolver.workflow_src.evolutionary_analysis.sequence_process import sequence_process
from enzyme_evolver.workflow_src.evolutionary_analysis.seqlogo import seqlogo
from enzyme_evolver.workflow_src.evolutionary_analysis.HeaderTrim import HeaderTrim
from enzyme_evolver.workflow_src.evolutionary_analysis.HeaderTrim import GetHeader
from enzyme_evolver.workflow_src.evolutionary_analysis.HeaderTrim import Make_Header
from enzyme_evolver.workflow_src.evolutionary_analysis.HeaderTrim import ReplaceHeader
from pathlib import Path
class Phylogeny:
    def __init__(self, folder_id, alignment, bootstrap_ancestral=False, tree_number = '50', test_mode=False, print_log=True):
        #masterpath is under /EnzymeEvolver
        self.master_path = Path.cwd()
        #path is under /EnzymeEvolver/enzyme_evolver/
        self.path = f"{Path(__file__).parents[1]}"
        self.folder_id = folder_id
        #folder_path is each jobs's folder
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.db_job = Job.objects(folder_id=folder_id)[0]
        self.test_mode = test_mode
        self.print_log = print_log
        self.alignment = alignment
        self.bootstrap_ancestral = bootstrap_ancestral
        self.tree_number = tree_number
        self._log(f"Folder ID = {self.folder_id}")
        self._log(f"Folder path = {self.folder_path}")

    def pre_process(self):
        self.db_job.update_status('pre-processing input sequence file', print_log=self.print_log)
        if self.alignment=='alignment':
            file = 'input_aln.fasta'
        else:
            file = 'input.fasta'
        print(file)
        sp.run(f'dos2unix {self.folder_path}/{file}', shell=True)
        spec_char_script = f"{self.path}/workflow_src/evolutionary_analysis/rm_specialChara.sh"
        sp.run(f'sh {spec_char_script} {self.folder_path}/{file}', shell=True)
        ##short the header
        HeaderTrim(f'{self.folder_path}/{file}')
        ##remove duplicate if it is seqeunce file
        if self.alignment=='sequence':
            Rm_duplicate = f'cd-hit -i {self.folder_path}/{file}  -o {self.folder_path}/{file}2 -c 1.0'
            sp.run(Rm_duplicate,shell=True)
            sp.run(f'mv {self.folder_path}/{file}2 {self.folder_path}/{file}',shell=True)
        ##get the original time
        GetHeader(f'{self.folder_path}/{file}', f'{self.folder_path}/original_header.txt')
        ## make unique header for original headers
        Make_Header(f'{self.folder_path}/original_header.txt', f'{self.folder_path}/uniq_header.txt')
        # replace original header with uniq header
        ReplaceHeader(f'{self.folder_path}/{file}', f'{self.folder_path}/original_header.txt',
                      f'{self.folder_path}/uniq_header.txt', f'{self.folder_path}/replace.txt')

        self.db_job.add_file('original_header.txt')
        self.db_job.add_file('uniq_header.txt')
        self.db_job.add_file('replace.txt')

    def muscle_alignment(self, sequence='input.fasta'):
        self.db_job.update_status('Multiple sequence alignment by muscle', print_log=self.print_log)
        ##seqeunce alignment
        sequence_process(ofolder=f"{self.folder_path}/", file=sequence, trim=False,
                         ofseq_ali='input_aln.fasta', ofhtml='input_aln_muscle.html')
        seqlogo(file='input_aln.fasta', ofolder=f"{self.folder_path}/")

        self.db_job.add_file('input.fasta')
        self.db_job.add_file('input_aln_muscle.html')
        self.db_job.add_file('input_aln_sequence_logo.png')


    def start_tree(self):
        self.db_job.update_status('starting tree generation', print_log=self.print_log)
        ####using  alignment to generate the starting tree (*.raxml.bestTree)

        sp.run(
            f'raxml-ng --msa {self.folder_path}/input_aln.fasta --prefix {self.folder_path}/start --model LG  '
            f'--seed 2 --tree pars{"{"}10{"}"},rand{"{"}10{"}"} --redo --threads 4 --force', shell=True)
        self.db_job.add_file(f'start.raxml.bestTree')
        self.db_job.add_file(f'start.raxml.mlTrees')



    def model_sele(self):
        self.db_job.update_status('selecting the best substitution model', print_log=self.print_log)
        # check the reasonability of the starting tree and evaulate the protein substitution model and output the top3_model.txt
        models = ['DAYHOFF','DCMUT','LG','JTT','MTREV','WAG','RTREV','CPREV','VT','BLOSUM62','MTMAM']
        for model in models:
            sp.run(f'raxml-ng --evaluate --msa {self.folder_path}/input_aln.fasta --model {model}+G+F '
                   f'--tree {self.folder_path}/start.raxml.bestTree --prefix {self.folder_path}/E_{model} '
                   f' --data-type AA --redo --force', shell=True)
        os.chdir(f'{self.folder_path}')
        sp.run(
            f'sh {self.master_path}/enzyme_evolver/workflow_src/evolutionary_analysis/model_sele.sh ./',
            shell=True)
        os.chdir(self.master_path)
        self.db_job.add_file(f'top3_model.txt')




    def final_tree(self):
        self.db_job.update_status('calculating the phylogenetic tree', print_log=self.print_log)
        with open(f'{self.folder_path}/top3_model.txt') as f:
            best_model = f.readline().strip()
            sp.run(f'raxml-ng --msa {self.folder_path}/input_aln.fasta --model {best_model}+G+F --prefix {self.folder_path}/T_{best_model} '
                   f'--seed 2 --tree pars{"{"}{self.tree_number}{"}"},rand{"{"}{self.tree_number}{"}"} --redo --force',shell=True)
        ReplaceHeader(f'{self.folder_path}/T_{best_model}.raxml.bestTree', f'{self.folder_path}/uniq_header.txt',
                      f'{self.folder_path}/original_header.txt', f'{self.folder_path}/replace.txt')
        self.db_job.update_status(f'the best substitution model is: {best_model}', print_log=self.print_log)
        self.db_job.add_file(f'T_{best_model}.raxml.bestTree')


    def final_tree_bs_ac(self):
        with open(f'{self.folder_path}/top3_model.txt') as f:
            best_model = f.readline().strip()
            #generate tree using best model
            sp.run(
                f'raxml-ng --msa {self.folder_path}/input_aln.fasta --model {best_model}+G+F --prefix {self.folder_path}/T_{best_model} '
                f'--seed 2 --tree pars{"{"}50{"}"},rand{"{"}50{"}"} --redo --force', shell=True)
            ##bootstrap
            self.db_job.update_status('bootstraping', print_log=self.print_log)
            sp.run(f'raxml-ng --bootstrap --bs-trees 100 --msa {self.folder_path}/input_aln.fasta --model {best_model}+G+F '
                   f' --prefix {self.folder_path}/bs_{best_model} --seed 2 --tree pars{"{"}20{"}"},rand{"{"}20{"}"} --redo --threads 4 --force',shell=True)
            ###supported values
            self.db_job.update_status('mapping suppoted values', print_log=self.print_log)
            sp.run(f'raxml-ng --support --tree {self.folder_path}/T_{best_model}.raxml.bestTree --bs-trees {self.folder_path}/bs_{best_model}.raxml.bootstraps'
                   f' --prefix {self.folder_path}/S_{best_model} --redo --threads 4 --force',shell=True)
            # ancestral sequence calculation
            self.db_job.update_status('calculating the ancestral seqeunce', print_log=self.print_log)
            sp.run(
                f'raxml-ng --msa {self.folder_path}/input_aln.fasta --model {best_model}+G+F --tree {self.folder_path}/T_{best_model}.raxml.bestTree '
                f'--ancestral --prefix {self.folder_path}/anc_{best_model} --data-type AA --redo --force --threads 4', shell=True)
        ReplaceHeader(f'{self.folder_path}/T_{best_model}.raxml.bestTree', f'{self.folder_path}/uniq_header.txt',
                      f'{self.folder_path}/original_header.txt', f'{self.folder_path}/replace.txt')
        ReplaceHeader(f'{self.folder_path}/S_{best_model}.raxml.support', f'{self.folder_path}/uniq_header.txt',
                      f'{self.folder_path}/original_header.txt', f'{self.folder_path}/replace.txt')
        ReplaceHeader(f'{self.folder_path}/anc_{best_model}.raxml.ancestralTree', f'{self.folder_path}/uniq_header.txt',
                      f'{self.folder_path}/original_header.txt', f'{self.folder_path}/replace.txt')
        self.db_job.add_file(f'T_{best_model}.raxml.bestTree')
        self.db_job.add_file(f'S_{best_model}.raxml.support')
        self.db_job.add_file(f'anc_{best_model}.raxml.ancestralProbs')
        self.db_job.add_file(f'anc_{best_model}.raxml.ancestralStates')
        self.db_job.add_file(f'anc_{best_model}.raxml.log')
        self.db_job.add_file(f'anc_{best_model}.raxml.ancestralTree')

    def clean_dir(self):
        self.db_job.update_status('cleaning the directory', print_log=self.print_log)
        sp.run(f'sh {self.master_path}/enzyme_evolver/workflow_src/evolutionary_analysis/clean.sh {self.folder_path}/',shell=True)

    def seq2tree(self):
        self.db_job.update_status('phylogenetic tree is being calculated', print_log=self.print_log)
        ##if input is only sequences without alignment, do alignment first
        seq_file = Path(f"{self.folder_path}/input.fasta")
        self.pre_process()
        if seq_file.is_file():
            self.muscle_alignment()
        self.db_job.add_file('input_aln.fasta')
        #remove special character in input_aln.fasta
        spec_char_script = f"{self.path}/workflow_src/evolutionary_analysis/rm_specialChara.sh"
        sp.run(f'sh {spec_char_script} {self.folder_path}/input_aln.fasta', shell=True)
        #alignment format: fasta
        self.start_tree()
        self.model_sele()
        if self.bootstrap_ancestral == False:
            self.final_tree()

            # sp.run(
            #     f'sh {self.master_path}/enzyme_evolver/workflow_src/evolutionary_analysis/phylo.sh {self.folder_path}/{self.alignment}  {self.folder_path}/',shell=True)




            # [batch_replace.batch_replace(ofolder + treefile, ofolder + 'replace.txt') for treefile in os.listdir(ofolder) if
            #  fnmatch.fnmatch(treefile, '*bestTree')]
        else:
            self.final_tree_bs_ac()
            # [batch_replace.batch_replace(ofolder + treefile, ofolder + 'replace.txt') for treefile in
            #  os.listdir(ofolder) if
            #  fnmatch.fnmatch(treefile, '*support')]
        # self.clean_dir()
        ###restore the header of input.fasta and input_aln.fasta. replace the uniq_header.txt with original_header.txt
        if self.alignment=='sequence':
            ReplaceHeader(f'{self.folder_path}/input.fasta', f'{self.folder_path}/uniq_header.txt',
                      f'{self.folder_path}/original_header.txt', f'{self.folder_path}/replace.txt')
        ReplaceHeader(f'{self.folder_path}/input_aln.fasta', f'{self.folder_path}/uniq_header.txt',
                      f'{self.folder_path}/original_header.txt', f'{self.folder_path}/replace.txt')


    def make_zip(self):
        self.db_job.update_status('Making zip file', print_log=self.print_log)
        os.chdir(f'{self.folder_path}')
        zipfile = self.folder_id + '.zip'
        sp.run(f'zip -r {zipfile} *', shell=True)
        os.chdir(self.master_path)
        self.db_job.add_file(f'{self.folder_id}.zip')

    def job_finished(self):
        self.db_job.update_status('Job finished', print_log=self.print_log)

    def _log(self, msg):
        if self.print_log is True:
            print(msg)
    def run(self):
        self.seq2tree()
        self.make_zip()
        self.job_finished()


if __name__ == '__main__':
    # from pathlib import Path
    # working_dir = str(Path(__file__).parents[2])
    #fileid = 'c4d88fa4-c2a4-4123-8997-46064890acc2'
    from enzyme_evolver.mongo.default_connection import make_default_connection
    make_default_connection()

    from enzyme_evolver.app.workflow.functions import job_functions
    folder_id = job_functions.create_new_job('test_protein', 'Phylogenetic Tree', no_user=True)

    # test_mode requires a blast.xml file in /database
    evol_analysis = EvolAnalysis(folder_id, test_mode=True, print_log=True)
    evol_analysis.save_fasta('test_protein', 'MG')
    evol_analysis.run()
    fileid = 'POR'
    alignment = 'enzyme_evolver/database/' + fileid + '/' + 'homologous_sequences90_mafft_trimmed.fas'
    ofolder = 'enzyme_evolver/database/' + fileid + '/'

    seq2tree(alignment, fileid, ofolder,  bootstrap_ancestral=True)
