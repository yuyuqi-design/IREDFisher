import subprocess as sp
import os
import sys
import fnmatch
import shutil
from enzyme_evolver.workflow_src.homology_modelling import auto_modeller
from enzyme_evolver.workflow_src.homology_modelling.modelling import model_class
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
from pathlib import Path

class HomoModelling():

    def __init__(self, folder_id, auto_single, auto_multiple, template_single, homodimer,start_residue='', end_residue='', test_mode=False, print_log=True):
        self.masterpath = Path.cwd()
        self.path = f"{Path(__file__).parents[1]}"
        self.folder_id = folder_id
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.db = f"{self.path}/database/pdb_95.pir"
        self.codes = f'{self.folder_path}/codes.txt'
        self.auto_single = auto_single
        self.auto_multiple = auto_multiple
        self.template_single = template_single
        self.homodimer = homodimer
        self.start_residue = start_residue
        self.end_residue = end_residue
        self.db_job = Job.objects(folder_id=folder_id)[0]
        self.db_job.update_notes(f"default = {auto_single}, use multiple template = {auto_multiple}, use your own template = {template_single}, "
                                 f"homodimer modelling = {homodimer}")
        self.test_mode = test_mode
        self.print_log = print_log

        self._log(f"Folder ID = {self.folder_id}")
        self._log(f"Folder path = {self.folder_path}")

    def create_sequence_files(self, fasta_string):
        self.save_fasta(fasta_string)
        self.mk_dir()

    def save_fasta(self, fasta_string):  # retrieve sequences from webpage input and save them as protein.fasta file

        fasta_filename = 'allseq.fasta'

        path = f"{self.folder_path}/{fasta_filename}"

        with open(path, 'w') as file:
            file.write(fasta_string)
            file.write('\n')
        self.db_job.add_file(fasta_filename)
        print(self.db)
        sp.run(f'dos2unix {path}', shell=True)

        spec_char_script = f"{self.path}/workflow_src/evolutionary_analysis/rm_specialChara.sh"
        sp.run(f'sh {spec_char_script} {self.folder_path}/allseq.fasta', shell=True)

        # split the fasta sequences into individual files and save all fasta titles into codes.txt file
        splitSeq_script = f"{self.path}/workflow_src/homology_modelling/split_sequences2.sh"
        # command_rm = f'dos2unix {self.folder_path}/{fasta_filename}'
        # print(command_rm)

        command_split = f'sh {splitSeq_script} {self.folder_path}/{fasta_filename} {self.folder_path}/ > {self.codes}'
        print(command_split)
        sp.run(f'sh {splitSeq_script} {self.folder_path}/{fasta_filename} {self.folder_path}/ > {self.folder_path}/codes0.txt', shell=True)
        sp.run(f"sort {self.folder_path}/codes0.txt |uniq > {self.codes}",shell=True)

    def mk_dir(self):  # create directory for each sequence
        with open(self.codes) as fin:
            codes = fin.readlines()
            for code in codes:
                code = code.strip()
                try:
                    # make subdirectory based on code
                    os.mkdir(f"{self.folder_path}/{code}")
                except FileExistsError:
                    shutil.rmtree(f"{self.folder_path}/{code}")
                    os.mkdir(f"{self.folder_path}/{code}")
                sp.run(f"cp {self.folder_path}/{code}.fasta {self.folder_path}/{code}", shell=True)

    def run(self, IREDFisher=False):
        with open(self.codes) as fin:
            codes = fin.readlines()
            # IREDFisher_template_file = open(f"{self.folder_path}/template.txt", "w")
            for i, code in enumerate(codes):
                self.db_job.update_status(f'Creating models for {i + 1}th of {len(codes)}', print_log=self.print_log)
                code = code.strip()
                subfolder_path = f"{self.folder_path}/{code}"
                os.chdir(subfolder_path)
                auto_modeller.fasta2pir(code)
                ######inspect the sequences and build models
                if IREDFisher:
                    auto_modeller.template_search(code, self.db)
                    templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain, template2_PdbCode, template3_chain, template3_PdbCode = auto_modeller.check_template(
                        code)
                    seqTemplate = [template1_PdbCode,template2_PdbCode, template3_PdbCode]
                    print(seqTemplate)
                    PDBIRED = ['3zgy', '4d3d', '4d3s', '4oqy', '4oqz', '5a9t', '5ocm', '5ojl', '6eod', '6jit', '6jiz','6grl','5g6r']
                    print(PDBIRED)
                    if any(x in seqTemplate for x in PDBIRED):
                        sp.run(f"cat {code}.fasta >> ../allseq_ired.txt", shell=True)
                        # IREDFisher_template_file.write(code + '    ' + template1_PdbCode +'\n')
                        sp.run(f"cat info.txt >> ../template.txt", shell=True)
                        self.create_models(code)
                    else:
                        sp.run(f"sed -i '/{code}$/d' ../codes.txt",shell=True)
                ##### no inspection of sequences
                else:
                    # pass
                    self.create_models(code)
                try:
                    #copy the model structure to folder ID
                    print(code)
                    sp.run(f'cp {code}.B99990001.pdb ../{code}.pdb', shell=True)
                    self.db_job.add_file(f'{code}.pdb')
                    sp.run(f'cat info.txt >> ../Modelling_template_summary.txt', shell=True)
                    # self.summary_modelling()
                    sp.run(f'sh {self.masterpath}/enzyme_evolver/workflow_src/homology_modelling/clean.sh', shell=True)
                    os.chdir(self.masterpath)
                    sp.run(f"rm -r {subfolder_path}", shell=True)
                except:
                    pass
            # IREDFisher_template_file.close()
            self.db_job.add_file(f'template.txt')
            self.db_job.add_file(f'Modelling_template_summary.txt')
            self.db_job.update_status(f'All models are finished', print_log=self.print_log)
        #sp.run(f'sh {self.path}/workflow_src/homology_modelling/best.sh {self.folder_path}/', shell=True)
        #sp.run(f'sh {self.path}/workflow_src/homology_modelling/rename.sh {self.folder_path}/', shell=True)
        #self.add_models()
        self.make_zip()
        self.job_finished()

    def create_models(self, seqcode):
        if self.auto_single == True:
            self.auto_single_modelling(seqcode)
        if self.auto_multiple == True:
            self.auto_multiple_modelling(seqcode)
        if self.template_single != False:
            sp.run(f'cp ../{self.template_single} ./', shell=True)
            self.template_single_modelling(seqcode)
        if self.homodimer == True:
            sp.run(f'cp ../dimerT.pdb ../singleT.pdb ./', shell=True)
            self.homodimer_modelling(seqcode)

    def auto_single_modelling(self,seqcode):
        #print(self.db)
        auto_modeller.template_search(seqcode, self.db)
        templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain, template2_PdbCode, template3_chain, template3_PdbCode = auto_modeller.check_template(
            seqcode)
        if templates:
            pass
            # auto_modeller.align2d(seqcode, template1, template1_PdbCode, template1_chain)
            # auto_modeller.trim_TerGap(str(seqcode).strip() + '-' + template1 + '.ali')
            # auto_modeller.build_model(seqcode, template1)
        else:
            pass

    def auto_multiple_modelling(self,seqcode):
        auto_modeller.template_search(seqcode, self.db)
        try:
            templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain, template2_PdbCode, template3_chain, template3_PdbCode = auto_modeller.check_template(
            seqcode)
            if templates and template1 and template2 and template3 :
                auto_modeller.salign(seqcode, template1_PdbCode, template1_chain, template2_PdbCode, template2_chain,
                                 template3_PdbCode, template3_chain)
                auto_modeller.align2d_mult(seqcode)
                # auto_modeller.trim_TerGap(str(seqcode).strip() + '-' + 'multi' + '.ali')
                auto_modeller.build_multi_models(seqcode, template1, template2, template3)
            else:
                pass
        except:
            pass

    def template_single_modelling(self,seqcode):
        template1 = self.template_single.split('.')[0] + 'A'
        template1_PdbCode = self.template_single.split('.')[0]
        template1_chain = 'A'
        auto_modeller.align2d(seqcode, template1, template1_PdbCode, template1_chain)
        auto_modeller.trim_TerGap(str(seqcode).strip() + '-' + template1 + '.ali')
        auto_modeller.build_model(seqcode, template1)

    def homodimer_modelling(self,seqcode):
        auto_modeller.single2dimer_ali(seqcode, self.start_residue, self.end_residue, '2')
        auto_modeller.build_chains(seqcode)

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
    folder_id = job_functions.create_new_job('test_protein', 'Homo Modelling', no_user=True)
    #folder_id = 'test'

    # fasta_string = '>Test_seq \nMDFATGFADKVALVVGGGRGIGSAVVEELARRGARVVVADTDTL' \
    #                'PSQYNHYQSTQVSGYADAQKLAARLTEEGLQVTAAQADATDEDQVSRLYADLAEQAGRLDVVVNAFGVTHVC' \
    #                'PVERMELAEFQRVVSGNLDGVFLSSKHAVPLLRDSGGGAIINFSSVSGRSGFAKVAHYCAGKFGVVGFTAALAQ' \
    #                'EVARDGIRVNAVCPGIVRSNMWRYLLSEFVRPGETEDECWERMRSMIPQREFQTPKDLAELVVYLAGATKVTGQAI' \
    #                'SVDGGMTAP\n>Test2_seq \nMRNMDFATGFADKVALVVGGGRGIGSAVVEELARRGARVVVADTDTLPSQYNHYQ' \
    #                'STQVSGYADAQKLAARFTEEGLRVTAAQADATDEDQVSRLYADLAEQAGRLDVVVNAFGVTHVCPVEKMELAEFQRVVS' \
    #                'GNLDGVFLSSKHAVPLLRDSGGGAIINFSSVSGRSGFAKVAHYCAGKFGVVGFTAALAQEVARDGIRVNAVCPGIVRSN' \
    #                'MWRYLLSEFVRPGETEDECWERMRSMIPQREFQTPKDLAELVVYLAGATKVTGQAISVDGGMTAP\n'

    fasta_string = '>Test2 \nMSEVAVIGLGRMGSALAKALITSGRSVTVWNRTPGKAEALEHLGASRAETPSAAIAASST' \
                   'LIVCLSDYAATSMVLDDACATDLLQGKTVVQLTSGTPKQARQLEEWVANRGGSYLDGAIS' \
                   'AWPSQIGGPEASIVIAGRESVFTSLEASLRLLAPNLTHVGNDISRAKVLFNAALAYFAGH' \
                   'WIGFSHGAAMCAAEGMDVAEFGETIASLSPMFADDLRHMGRAIEGNRFADPQSTIRSVGV' \
                   'DITRLVEIADDLNINTAFPAFASDLFRNATDAGYGAEEHCAIVKVIRAW'
    hm = HomoModelling(folder_id)
    hm.create_sequence_files(fasta_string)
    hm.run()
