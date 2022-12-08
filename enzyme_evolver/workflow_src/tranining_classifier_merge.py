import pandas as pd
def classifier_concatenate(df_classifier,training_x = 'training_noclassifier.csv',training_xy = 'training_classifier.csv',*columns):
    df = df_classifier
    df_trainingx = pd.read_csv(training_x)
    concatenate_classifier = []
    for column in columns:
        print(column)
        classifier = df[column].values.tolist()
       # print(len(classifier))
        concatenate_classifier= concatenate_classifier + classifier
    df_trainingx['classifier'] =pd.Series(concatenate_classifier)
    df_trainingxy = df_trainingx
    df_trainingx.to_csv(training_xy, index=False)
    return df_trainingxy





if __name__ == '__main__':
    folder = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter'
    input_file = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/cutoff0.5/IRED_SLM_classifier.csv'
    training_x = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/training_noclassifier.csv'
    training_xy = '/hdd/data/validation_case/paper_sciAdv/result2/hyperparameter/cutoff0.5/training_classifier.csv'
    classifier_concatenate(input_file,training_x,training_xy, 'cpd3_classifier','cpd4_classifier','cpd5_classifier','cpd6_classifier','cpd10_classifier','cpd11_classifier')