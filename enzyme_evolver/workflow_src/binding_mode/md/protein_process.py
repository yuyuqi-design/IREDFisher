import subprocess as sp
def protein_process(file):
    # protonation = 'pdb2pqr --ff=AMBER --ffout=AMBER --ph-calc-method=propka --with-ph=7.4  ' + file  + ' ' + file +'.pqr ' + '--summary'
    # sp.run(protonation,shell=True)
    with open(file + '.pqr') as pqr:
        name = file.split('.')[0]
        with open( name + '_propka.pdb', 'w') as pdb:
            content = pqr.readlines()
            
            pdb.writelines(content[:len(content)-1])


protein_process('P95480_best_model.pdbqt')