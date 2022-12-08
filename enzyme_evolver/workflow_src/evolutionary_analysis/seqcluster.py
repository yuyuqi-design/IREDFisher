###cluster the fasta sequences by threshold 90%
#require cd-hit installed

import subprocess as sp
import sys

def seqcluster(seq='homologous_sequences.fas'):
    seq = str(seq).strip()
    #print(seq)
    command_cluster = 'cd-hit -i '+ seq +' -o ' + seq.split('.')[0] + ' -T 2 -c 0.90'
    sp.run(command_cluster, shell=True)
    #command_getRep = "grep '\*' " + seq.split('.')[0]+".clstr|awk -F \> '{print $2}'|awk -F . '{print $1}'|sort > seq_rep_id_" + seq.split('.')[0]+".txt"
    #sp.run(command_getRep, shell=True)
    #command_catsequence = "sh catsequence90.sh " + "seq_rep_id_" + seq.split('.')[0]+".txt"
    #sp.run(command_catsequence, shell=True)
    command_cluster_organism = "grep \> homologous_sequences|awk -F \> '{print $2}' > cluster_organism.txt"
    sp.run(command_cluster_organism,shell=True)
    sp.run('rm *fasta', shell = True)
    
    
    # with open('homologous_clustered90.fasta') as output:
    #     with open("seq_rep_id_" + seq.split('.')[0]+".txt") as clustered_seq:
    #         Ids = clustered_seq.readlines()
    #         for id in Ids:
