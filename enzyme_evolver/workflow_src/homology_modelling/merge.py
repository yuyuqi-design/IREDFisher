import pandas as pd
def merge(csv1, csv2):
    df_rmsd = pd.read_csv(csv1)
    df_modelling = pd.read_csv(csv2,delimiter='\t',names=['Sequence', 'Template_code', 'SequenceIdentity','Model','Function','Organism'])
    df = pd.merge(df_modelling,df_rmsd, on='Model')
    return df

if __name__ == "__main__":
    csv1 = 'models_rmsd2ref.csv'
    csv2 = 'Modelling_template_summary.txt'
    colums = ['Sequence', 'Template_code', 'SequenceIdentity','Model','Function']
    df = merge(csv1, csv2)
    df.to_csv('test.csv',index=False)


