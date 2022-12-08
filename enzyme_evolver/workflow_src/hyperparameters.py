#file_best_binding mode, Vina score, Acidic residue number, Basic residue number, His apprearence,classifier
#AspRedAm_cpd_3_best_mode_009,  -4.2,    0,   0,   0,   1
import pandas as pd
import numpy as np
# def tune_hyperparameters(file_name, vinaScore,acidic,basic, his,a4vinascore, b4acidic, c4basic, d4his):
#     refined_score = a4vinascore * vinaScore + b4acidic * acidic + c4basic * basic + d4his * his
#     return file_name, refined_score, a4vinascore, b4acidic, c4basic, d4his

def run_tuning(df_trainingxy,accurary_csv):
    accuracy_parameters = []

    for a4vinascore in np.arange(2,5,1):
        a4vinascore = "{:.2f}".format(a4vinascore)
        for b4acidic in [1]:
            b4acidic = "{:.2f}".format(b4acidic)
            for c4basic in np.arange(3,10,1):
                c4basic = "{:.2f}".format(c4basic)
                for d4his in np.arange(-10,-8,1):
                    df = df_trainingxy
                    d4his = "{:.2f}".format(d4his)
                    # print('type')
                    # print(type(a4vinascore))
                    # print(b4acidic)
                    # print(c4basic)
                    # print(d4his)
                    # print(df['Vina_score'] )
                    # print(df['Acidic_residue_number'])
                    # print(df['Basic_residue_number'])
                    # print(df['His_apprearence'])
                    df['refined_score'] = float(a4vinascore) *  df['Vina_score'] + float(b4acidic) * df['Acidic_residue_number'] + float(c4basic) * df['Basic_residue_number']  + float(d4his) * df['His_apprearence']
                    #print(df['refined_score'])
                    df['a4vinascore'] = a4vinascore
                    df['b4vinascore'] = b4acidic
                    df['c4basic'] = c4basic
                    df['d4his'] = d4his
                    df = df.sort_values(by=['refined_score'])


                    #total positives
                    N_positive_total = df['classifier'].sum()
                    ## top20 positives
                    N_positive_top20 = df[:20]['classifier'].sum()
                    ############whole ranking####################
                    # Total = len(df)
                    # picked_True_positive = df[:N_positive]['classifier'].sum()
                    # accurary = float(picked_True_positive) / float(N_positive)
                    #refined_scores_file = str(accurary) + '_' +str(a4vinascore) + '_' + str(b4acidic) +'_' + str(c4basic) +'_' + str(d4his) + 'refined_score.csv'
                    ############whole ranking####################
                    accurary = float(N_positive_top20)/float(N_positive_total)
                    accurary_parameter = [accurary,a4vinascore, b4acidic,c4basic,d4his]
                    print(accurary_parameter)
                    #df.to_csv(refined_scores_file, index=False)
                    accuracy_parameters.append(accurary_parameter)
                    
    results = pd.DataFrame(accuracy_parameters, columns=['accurary','a4vinascore','b4acidic','c4basic','d4his']).sort_values(by=['accurary'])
    #print(results)
    results.to_csv(accurary_csv, index=False)



if __name__ == '__main__':
    folder = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter'
    input_file = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/training.csv'
    #refined_scores_file = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/refined_scores.csv'
    output_file = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/accuracy.csv'
    run_tuning(input_file,output_file)



