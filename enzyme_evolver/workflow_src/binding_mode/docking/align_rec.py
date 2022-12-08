#####################################################################################
def align(receptors_file):
#align all the proteins in receptors_file
    with open(receptors_file, 'r') as recs_file:
        all_proteins = recs_file.readlines()
        refer_pro = all_proteins[0].strip().split('.')[0]
        #print(all_proteins)
        for protein in all_proteins:
            cmd.load(protein.strip())
        cmd.alignto(refer_pro)
        for protein in all_proteins:
            name = protein.strip().split('.')[0]
            cmd.save('aln_' + name + '.pdb', name)