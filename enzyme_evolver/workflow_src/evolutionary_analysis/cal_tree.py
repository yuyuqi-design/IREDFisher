#run phylip tree generation by maximum likelihood algorithm
import pandas as pd 
from batch_replace import batch_replace
# def cal_tree(infile,input):
#     runtree = 'phylip proml <' + input
#     sp.run(runtree,shell=True)

def node_name(treefile,namefile, seqfile):
    #namefile: seq_tax_refine.csv: ID,Organism,Tax_kingdom,Tax_group,Unnamed: 0,Accession,Description,Sequence,Tax_blastName
    df = pd.read_csv(namefile)
    dfseq = pd.read_csv(seqfile)
    dfseq.sort_values(by='sequence_name', inplace=True)
    df['sequence_name'] = df['Sequence'].str.split('\n',expand=True)[0].str.strip().str.split('>',expand=True)[1]
    df['node_name'] = df['Organism'] + '      ' + df['Tax_blastName'].str.strip()
    df.sort_values(by='sequence_name', inplace=True)
    dfnew = pd.merge(df,dfseq, on='sequence_name')
    ##replace seqnumber in tree file with node_name (organism     family)
    dfnew[['seqnumber','node_name']].to_csv('replace.txt',sep=';',index=False,header=None)
    batch_replace(treefile,'replace.txt')

    #replace seqnumber in tree file with organism
    # replace_list = list(dfnew['seqnumber'] + ';' +  dfnew['Organism']
    

    # f = open(treefile)
    # tree = f.read()
    # for line in replace_list:
    #     stxt = line.split(';')[0]
    #     # print(stxt)
    #     rtxt = line.split(';')[1]
    #     # print(rtxt)
    #     tree = tree.replace(stxt.strip() + ':', rtxt.strip()+':')
    # with open('outtree_taxno.txt', 'w') as output:
    #     output.write(tree)
    # f.close()

treefile = 'outtree'  
namefile = 'seq_tax_refine.csv'
seqfile = 'temp2.phy'
node_name(treefile,namefile, seqfile)