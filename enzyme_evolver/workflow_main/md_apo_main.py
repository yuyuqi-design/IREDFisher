from enzyme_evolver.workflow_src.binding_mode.docking.mkdir import mk_dir
from enzyme_evolver.workflow_src.binding_mode.md import prepare_md
import sys
import subprocess as sp
def md_apo():
    sp.run('cp enzyme_evolver/workflow_src/binding_mode/md/*mdp ' + folder, shell=True)
    sp.run('cp enzyme_evolver/workflow_src/binding_mode/md/job0.sh ' + folder,shell=True)
    name_dir = mk_dir(rec_file,lig_file,folder)
    print(name_dir)
    sub_folder = folder + name_dir + '/'
    sp.run(
        'cp ' + folder + lig_file + ' ' + folder +  rec_file + ' ' + folder + name_dir, shell=True)
