import subprocess as sp
def complex_process(pro_file,lig_file):
    with open(pro_file) as pro:
        with open(lig_file) as lig:
            with open('complex.pdb', 'w') as output:
                pro_content = pro.readlines()
                lig_content = lig.readlines()
                complex_content = pro_content + lig_content
                output.writelines(complex_content)
    tleap = 'tleap -f tleap.in'
    sp.run(tleap, shell=True)
complex_process('P95480_best_model_propka.pdb', 'lig_out.pdb')
        


