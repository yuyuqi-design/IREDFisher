from math import sqrt
import subprocess as sp
import pandas as pd
import numpy
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver import database_functions
from pymol import cmd,stored
from pathlib import Path

class IREDFisher_scoring():
    def __init__(self, folder_id, N_atom_string='N01 UNK'):
        self.folder_id = folder_id
        self.N_atom_string=N_atom_string
        self.path = f"{Path(__file__).parents[1]}"
        self.folder_path = f"{self.path}/database/{folder_id}"
        self.masterpath = Path.cwd()
        # self.db_job = Job.objects(folder_id=folder_id)[0]

##

##calculate the distance of two points
    def distance_two_atoms(self, p1,p2):
    #p1,p2  both are list
        return sqrt((float(p1[0])-float(p2[0]))**2 + (float(p1[1])-float(p2[1]))**2 + (float(p1[2])-float(p2[2]))**2)
    ##Pymol: get residues around 8 angstrom of ligname of a pdbfile
    def get_residues(self,path, file='IR-29_1a_best_mode.pdb',ligname='UNK',best_number=2):
        stored.list=[]
        #load the structure file
        cmd.load(f"{path}/{file}")
        #object name of the structure file in pymol
        if len(cmd.get_object_list())==1:
            obj = file.split('.')[0]
        elif int(best_number) <10:
            obj = file.split('.')[0] + '_000'+str(best_number)
        else:
            obj = file.split('.')[0] + '_00'+str(best_number)
        #select ligand from the pdb file by residue name ligname
        cmd.select('lig', obj +' and resn '+ligname)
        #select r for residues within 8 angstrom of ligand
        cmd.select('r', 'lig around 8')
        #select only named alpha atoms from r and save them into stored.list
        cmd.iterate("(name ca and r)","stored.list.append((resi,resn))")
        #count basic and acidic and histidine residue numbers
        basic_number = 0
        acid_number = 0
        his_number = 0
        tyr_number = 0
        for i in stored.list:
            #print(i[1])
            #if i[1] == 'ASP' or i[1] == 'GLU' or i[1] == 'HIS' or i[1] == 'TYR':
            #if there is aspartate or glutamate,acidic number is added by 1
            if i[1] == 'ASP' or i[1] == 'GLU':
                acid_number = acid_number + 1
            # if there is lysine or arginine, basic number is added by 1
            if i[1] == 'LYS' or i[1] == 'ARG':
                basic_number = basic_number + 1
            #if there is a histidine, histidine number is added by 1
            if i[1] == 'HIS':
                his_number = his_number + 1
            if i[1] == 'TYR':
                tyr_number = tyr_number + 1
        #delete all objects from pymol
        cmd.delete('all')
        #return all acidic, basic, histidine residue number and the list of all residues
        return acid_number, basic_number, his_number, tyr_number, stored.list

    ##calculate distances between NADP C4 atom and every docking pose's atom for reaction (N 01) and return the first distance, the minimum distance + its corresponding rank, and all distance values for every pose.
    #input the docking output pdb file
    #string1 is for NAP C4 atom name string
    #string2 is for active atom for the ligand
    # file is the input pdb composed of protein, NAP and every docking poses
    def distance_active_atoms(self, string1 = 'C4N NAP',file='IR-92_1a_best_mode.pdb'):
        #get the coordinate of C4 atom of NAP from the structure file
        get_p1_line = f"grep --no-filename {string1} {file}"
        p1_line = sp.run(get_p1_line, capture_output=True, encoding="utf-8",shell=True)
        p1 = p1_line.stdout[32:54].split()
        print(p1)
        #get every active atom's coordinate in each docking pose
        get_p2_line = f"grep --no-filename {self.N_atom_string} {file}"
        p2_line = sp.run(get_p2_line, capture_output=True, encoding="utf-8",shell=True)
        p2s = p2_line.stdout
        print(p2s.split('\n')[0])
        #make a list to save every distance between C4 atom of NAP and the active atom in each docking pose.
        all_dist = []
        for line in p2s.split('\n'):
            if line:
                p2 = line[32:54].split()
                all_dist.append(self.distance_two_atoms(p1,p2))
        #turn every element in the list into float type
        all_dist_float = [float(i) for i in all_dist]
        print(all_dist)
        if all_dist:
        #get the distance of the first docking pose
            first = all_dist[0]
        #get the minimum distance between C4 atom of NAP and the active atom in each docking pose.
            minimum = min(all_dist_float)
        #get the rank of the minimum distance
            index_min = all_dist_float.index(min(all_dist_float))
        #return the first distance, minimum distance, the rank of the minimum disntance and distances for all docking poses
            return first,minimum, index_min, all_dist_float
        else:
            first = '999999'
            minimum = 999999
            index_min = 999999
            return first, minimum, index_min, all_dist_float



    ###get docking score of each pose from a pdbqt file
    def vina_scores(self, file='IR-92_1a_best_mode.pdbqt'):
        #get the lines containing the score value
        get_score_line = f"grep --no-filename 'REMARK VINA RESULT' {self.folder_path}/{file}"
        score_line = sp.run(get_score_line, capture_output=True, encoding="utf-8",shell=True)
        score_lines = score_line.stdout
        #print(score_lines)
        #make a list to save scores
        all_score = []
        for line in score_lines.split('\n'):
            if line:
                score = line.split()[3]
                #print(score)
                all_score.append(score)
        #turn element in score list into float type
        all_score_float = [float(i) for i in all_score]
        return all_score_float
    def combine_score_files(self):
        with open(f'{self.folder_path}/ligands.txt') as ligs:
            ligands = ligs.readlines()
            first_lig = ligands[0].split(".")[0] + '_sort_rescore.csv'
            print(first_lig)
            df_index = pd.read_csv(f'{self.folder_path}/{first_lig}')
            for lig in ligands[1:]:
                name = lig.split(".")[0]
                score_file = name + '_sort_rescore.csv'
                df_score = pd.read_csv(f'{self.folder_path}/{score_file}')
                df_index = pd.merge(df_index, df_score, right_index=True, left_index=True)
        df_index.to_csv(f'{self.folder_path}/all_scores.csv')

    def run_rescore(self):
        with open(f'{self.folder_path}/receptors.txt') as fin:
            with open(f'{self.folder_path}/ligands.txt') as fin2:
                recs = fin.readlines()
                ligs = fin2.readlines()
                for lig in ligs:
                    list_summary = []
                    for rec in recs:
                        lig = lig.split('.')[0]
                        file = rec.split('.')[0] + f'_{lig}_best_mode.pdb'
                        if Path(self.folder_path +'/' +file).is_file():
                            first, minimum, index_min, all_dist_float = self.distance_active_atoms(string1='C4N NAP', file=self.folder_path +'/' +file)
                            print(first)
                            if all_dist_float:
                                print(all_dist_float)
                                all_score_float = self.vina_scores(file + 'qt')
                                data_tuples = list(zip(all_dist_float, all_score_float))
                                df = pd.DataFrame(data_tuples, columns=['Distance', 'Score'])
                                # print(df)
                                df_dist = df[(df['Distance'] >= 3.5) & (df['Distance'] <= 6.0)]
                                # print(df_dist)
                                df_dist_sort = df_dist.sort_values(by='Distance')
                                # print(df_dist_sort)
                                Distance = (df_dist_sort.iloc[:1]['Distance']).to_string()
                                # print(Distance)
                                if 'Series' in Distance:
                                    Distance = '0 66'

                                Score = (df_dist_sort.iloc[:1]['Score']).to_string()
                                print(Score)
                                if 'Series' in Score:
                                    Score = '0 10'
                                best_mode_number = str(int(Distance.split()[0]) + 1)
                                Acid_number, Basic_number, His_number, Tyr_number, Stored_list = self.get_residues(self.folder_path, file, 'UNK', best_mode_number)
                                if His_number > 0:
                                    His_number = 1
                                    His_apprearence = '1'
                                if His_number == 0:
                                    His_number = -5
                                    His_apprearence = '0'

                                # print(His_number)
                                # raw tuple : file_best_binding mode; Docking score; Distance; Acidic residue number, Basic residue number
                                # tuple_data = (file.split('.')[0]+'_00'+best_mode_number,float(Score.split()[1]),float(Distance.split()[1]),float(Acid_number),float(Basic_number), float(His_number))

                                ##refined docking score output
                                if len(Distance)==1:
                                    best_mode = file.split('.')[0]
                                    refine_score = 4.0 * float(Score.split()[1]) + 1.0 * float(Acid_number) - 9.0 * float(His_number) + 9.0 * float(Basic_number)
                                    refine_score = "%.2f" % refine_score
                                    picking_pose_score = rec.split('.')[0] + ',    ' + best_mode + ',   ' + str(refine_score) + '\n'
                                else:
                                    best_mode = file.split('.')[0] + '_00' + best_mode_number
                                    refine_score = 4.0 * float(Score.split()[1])+ 1.0 * float(Acid_number) - 9.0 * float(His_apprearence) + 9.0 * float(Basic_number)
                                    # refine_score =  float(Score.split()[1]) -1.5 * float(Acid_number) - 2 * float(
                                    #     His_apprearence) + float(Basic_number)
                                    refine_score = "%.2f" % refine_score
                                # print('vina score is ' + Score.split()[1])
                                # print(Acid_number)
                                # print(His_number)
                                # print(Basic_number)
                                    picking_pose_score = rec.split('.')[0] +',    ' + best_mode + ',   ' + str(refine_score)  + '\n'
                                    print(picking_pose_score)
                                #tuple_data = (file.split('.')[0] + '_00' + best_mode_number,
                                              # str(float(Score.split()[1]) - 1.5 * float(Acid_number) - 2 * float(His_number) + float(
                                              #     Basic_number)))
                                # print(file)
                                    print(Stored_list)
                            else:
                                picking_pose_score = rec.split('.')[0] +',    ' + 'NONE' + ',   ' + str(999999)  +'\n'
                            list_summary.append(picking_pose_score)

                        else:
                            pass
                ###list_summary: IRED, Best_mode, refined score
                    list_summary2 = [x.split(',') for x in list_summary]
                    df_list_summary2 = pd.DataFrame(list_summary2, columns=['IRED', 'Best mode', 'Refined score'])
                    df_list_summary2['Refined score'] = df_list_summary2['Refined score'].astype(float)
                    df_list_summary3 = df_list_summary2.sort_values(by='Refined score')
                    df_list_summary3.to_csv(f"{self.folder_path}/{lig}_sort_rescore.csv", index=False)
                    IRED_top25 = df_list_summary3['IRED'].tolist()[:26]
                    for IRED in IRED_top25:
                        sp.run(f"cat {self.folder_path}/{IRED}.fasta >> {self.folder_path}/{lig}_optimized_panel_sequence.fasta",shell=True)
        self.combine_score_files()

    # PyMOL not running, entering library mode (experimental)
    #('IR-01_7b_best_mode_009', -9.0)
                    # with open(f'{self.folder_path}/{lig}_rescore.csv', 'w') as fout:
                    #     # for each in list_summary:
                    #     #     each = ',    '.join(each)
                    #     #     fout.write(each+'\n')
                    #     fout.write('IRED,    Refined_score\n')
                    #     fout.writelines(list_summary)
                    # #[print(i) for i in list_summary]
                    # sort_rescore = f"head -1 {self.folder_path}/{lig}_rescore.csv > {self.folder_path}/{lig}_sort_rescore.csv; " \
                    #                f"sed 1d {self.folder_path}/{lig}_rescore.csv|sort -n -k2 >> {self.folder_path}/{lig}_sort_rescore.csv "
                    # sp.run(sort_rescore,shell=True)

                    # self.db_job.add_file(f'{lig}_sort_rescore.csv')
                    # self.db_job.add_file(f'{lig}_optimized_panel_sequence.fasta')



if __name__ == '__main__':

    folder_id = 'c7584321-e9c5-4516-85a0-b31dc9604d93'
    Rescore = IREDFisher_scoring(folder_id)
    Rescore.run_rescore()
