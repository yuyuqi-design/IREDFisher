###only works when the ligand startwith 'HETATM' in PDB file. only one ligand in pdb fil is recommended.
import subprocess as sp

def complex_split(complex='complex.pdb', pro='pro.pdb', lig='lig.pdb'):
    remove_ANISO = f"sed -i '/ANISOU/d;/CONECT/d' {complex}"
    print(remove_ANISO)
    sp.run(remove_ANISO, shell=True)
    #get ligand by matching the word HETATM
    extract_lig = f"grep HETATM {complex} > {lig}"
    print(extract_lig)
    sp.run(extract_lig, shell=True)
    # get ligand by not matching the word HETATM
    extract_pro =  f"grep -v HETATM {complex} > {pro}"
    sp.run(extract_pro, shell=True)