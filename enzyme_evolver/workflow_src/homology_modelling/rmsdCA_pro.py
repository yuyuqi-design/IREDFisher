from pymol import cmd
import subprocess as sp
import csv
import pandas as pd


def rmsd_super_CA(recs='IsaList',folder=''):
    #align all the proteins in receptors_file
    rmsd = ['0.0']
    refer_pro = recs[0].strip().split('.')[0]
    # print(recs)
    cmd.load(folder+'/'+'pro.pdb')
    cmd.select('proca', 'pro and name CA')

    for protein in recs[1:]:
        # print(folder+'/'+protein.strip())
        cmd.load(folder+'/'+protein.strip())
        protein_name = protein.strip().split('.')[0]
        cmd.select(protein_name+'ca', protein_name + ' and name CA')
        rmsd_protein = cmd.super(protein_name+'ca', 'proca', cycles=0)
        rmsd_pro = "{0:.1f}".format(rmsd_protein[0])
        rmsd.append(rmsd_pro)
    cmd.delete('all')
    #print(rmsd)
    #print(rmsd)


    with open(f'{folder}/models_rmsd2ref.csv', 'w', newline='') as rmsdfile:
        writer = csv.writer(rmsdfile)
        writer.writerow(['Model','RMSD'])
        recs = [line.strip() for line in recs]
        writer.writerows(zip(recs,rmsd))
    df_rmsd = pd.read_csv(f'{folder}/models_rmsd2ref.csv')
    df_rmsd = df_rmsd.sort_values(by=['RMSD'])
    df_rmsd.to_csv(f'{folder}/models_rmsd2ref.csv',index=False)
    # print(df_rmsd)
    return df_rmsd


def split_positive_negative(recs='IsaList',folder='',rmsd=2.0):
    df_rmsd = rmsd_super_CA(recs,folder)
    postive = df_rmsd[df_rmsd['RMSD'] <= rmsd]
    negative = df_rmsd[df_rmsd['RMSD'] > rmsd]
    postive.to_csv(f'{folder}/positive.txt', index = False, header=False)
    negative.to_csv(f'{folder}/negative.txt', index=False, header=False)
    with open(f'{folder}/positive.txt') as pos:
        seqs = pos.readlines()
        for seq in seqs:
            code = seq.split('.')[0]
            pos_seq_command = f'cat {folder}/{code}.fasta >> {folder}/positive.fasta'
            sp.run(pos_seq_command,shell=True)
    with open(f'{folder}/negative.txt') as neg:
        seqs = neg.readlines()
        for seq in seqs:
            code = seq.split('.')[0]
            neg_seq_command = f'cat {folder}/{code}.fasta >> {folder}/negative.fasta'
            sp.run(neg_seq_command,shell=True)
    return postive, negative

if __name__ == '__main__':
    folder = '/home/g02808yy/Downloads/test'
    with open('/home/g02808yy/Downloads/test/receptors.txt') as f:
        recs = f.readlines()
        rmsd = split_positive_negative(recs=recs,folder=folder)
        #print(rmsd)