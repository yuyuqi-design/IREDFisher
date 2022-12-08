import pandas as pd
def classifier_generation(activity_csv_file,cutoff=0.2):
    df_activity = pd.read_csv(activity_csv_file)
    columns = df_activity.columns
    enzyme = columns[0]
    new_columns = [enzyme]
    for column in columns[1:]:
        new_column = column +'_classifier'
        new_columns.append(new_column)
        df_activity[new_column] = df_activity[column].apply(lambda x: 1 if x >= cutoff else 0)
    return df_activity[new_columns]


if __name__ == '__main__':
    folder = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter'
    activity_csv_file = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/IRED_SLM.csv'
    #refined_scores_file = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/refined_scores.csv'
    classifier_csv_file = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/cutoff0.5/IRED_SLM_classifier.csv'
    classifier_generation(activity_csv_file,classifier_csv_file,cutoff=0.5)

