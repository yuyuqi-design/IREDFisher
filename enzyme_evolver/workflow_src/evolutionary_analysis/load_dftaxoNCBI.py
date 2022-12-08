import pandas as pd

def load_dftaxoNCBI(path='enzyme_evolver/database/fullnamelineage.dmp'):
    print('LOADING FULLLINEAGE.DMP...')
    df_taxoNCBI = pd.read_csv(path, delimiter='|', names = ['TaxID','Organism','Tax','nan'])
    print('FULLLINEAGE.DMP LOADED OK')
    # remove tab from organism and tax collum
    print('PARSING DF..')
    df_taxoNCBI['Organism'] = df_taxoNCBI['Organism'].str.strip()
    df_taxoNCBI['Tax_kingdom'] = df_taxoNCBI['Tax'].str.strip().str.split(';', expand=True)[1]
    df_taxoNCBI['Tax_group'] = df_taxoNCBI['Tax'].str.strip().str.split(';', expand=True)[4]
    print('DF PARSED')
    df_taxoNCBI.to_csv('taxoNCBI.csv')
load_dftaxoNCBI()
