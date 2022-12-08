#Copyright (c) 2020 Dr.Yuqi Yu, Nigel Scrutton Group, Manchester Institute of Biotechnology, The University of Manchester, UK
# email to: yuqi.yu@manchester.ac.uk or yuyuqihappy@gmail.com
####################################################################################################################################################################
#function: replace a set of strings with another set of strings       
#Input: your file: file.txt, a replace.txt with each line formatted  as apple;banana (replace apple with banana)                                                                                                                   
#Usage: python3 batch_replace.py file.txt replace.txt
#output: file.txt_replaced
#####################################################################################################################################################################

import pandas as pd
def replace_extract(df='',ofolder=''):
    df_replace = df['Organism'] + ';' + df['Organism']+ '    '+df['Tax_blastName']
    df_replace.to_csv(ofolder+'replace.txt',header=False, index=False)

def batch_replace(file='file.txt', s_r_list='replace.txt',ofolder=''):
    #s_rlist_file content: apple;banana (replace apple with banana) 
    f = open(file)
    #read the file as f_content (string)
    f_content = f.read()
    #read the replace.txt as dataframe with delimiter as ; and name the first column as stxt (searching string), the second is rtxt (replace string)
    df = pd.read_csv(s_r_list, sep=';',names=['stxt','rtxt'])
    df.head()
    #get the length of stxt then add them as a new column into dataframe
    df['length'] = df['stxt'].str.strip().str.len()
    #Rank the dataframe according to the length of stxt in decreasing order. This can avoid the substring was replaced if it is in ascending order.
    df.sort_values(by='length',ascending=False,inplace=True)
    #turn the stxt and rtxt back to list
    replace_list = list(df['stxt'] + ';' + df['rtxt'])
    #Do the replacement
    for line in replace_list:
        stxt = line.split(';')[0]
        # print(stxt)
        rtxt = line.split(';')[1]
        # print(rtxt)
        f_content = f_content.replace(stxt.strip() , rtxt.strip())
    #write to a new file
    with open(ofolder+file + '_taxo', 'w') as output:
        output.write(f_content)
    f.close()


if __name__ == '__main__':
    batch_replace(sys.argv[1],sys.argv[2])
    
