from enzyme_evolver.workflow_src.classifier import classifier_generation
from enzyme_evolver.workflow_src.tranining_classifier_merge import classifier_concatenate
from enzyme_evolver.workflow_src.hyperparameters import run_tuning
import numpy as np
import os

def run_tune_cutoff(activity_csv, training_noclassifier_csv,training_xy_csv,cutoff, accuracy_csv, *columns):
    #activity_csv: input file
    df_classifier = classifier_generation(activity_csv,cutoff=cutoff)
    #print(df_classifier)
    df_trainingxy = classifier_concatenate(df_classifier, training_noclassifier_csv, training_xy_csv, *columns)
    run_tuning(df_trainingxy,accuracy_csv)


if __name__ == '__main__':
    for cutoff in np.arange(0.5,1.0,0.5):
        cutoff = float("{:.2f}".format(cutoff))
        for cpd in ['13']:
        # for cpd in ['8','9','13','14']:
            folder = f'/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/cpd{cpd}'
            working_folder = f'/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/cpd{cpd}/cutoff{cutoff}'
            try:
                os.mkdir(working_folder)
            except FileExistsError:
                print('directory exists.')
            activity_csv = f'{folder}/IRED_SLM_noP19.csv'
            training_noclassifier_csv = f'{folder}/cpd{cpd}vinaScore_distance_Nresidues.csv'
            training_xy_csv = f'{working_folder}/training.csv'
            accuracy_csv = f'{working_folder}/accuracy.csv'
            #run_tune_cutoff(activity_csv,training_noclassifier_csv,training_xy_csv,cutoff,accuracy_csv,'cpd3_classifier','cpd4_classifier','cpd5_classifier','cpd6_classifier','cpd10_classifier','cpd11_classifier')
            run_tune_cutoff(activity_csv, training_noclassifier_csv, training_xy_csv, cutoff, accuracy_csv, f'cpd{cpd}_classifier')
