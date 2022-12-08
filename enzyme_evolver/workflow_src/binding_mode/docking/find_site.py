import subprocess as sp
import os
import fnmatch
import pymol
from pymol import cmd
def site_calculate(structure,folder):
    #input structure and directory, output traj.clust, traj_analysis.dat, test_26_patch.pdb(location blob), test_26_plas.pdb(residues)
    epos_commd= 'EPOS_x86 -file '+ folder + structure +' -cluster 75 1 ' + folder +'traj.clust -analyse ' + folder + 'traj_analysis.dat'
    print(epos_commd)
    sp.run(epos_commd,shell=True)

def box_center(reference,folder):
    # if user uses ligand to specify the location of pocket
    cmd.load(folder+reference)
    box_center = cmd.centerofmass(reference.split('.')[0])
    cmd.delete('all')
    return box_center
def centers_patch(folder):

    centers = [box_center(patch,folder) for patch in os.listdir(folder) if fnmatch.fnmatch(patch,'*patch.pdb')]
    return centers

if __name__ == '__main__':
    folder = 'enzyme_evolver/database/site/'
    structure = 'test.pdb'
    site_calculate(structure, folder)
    centers = centers_patch(folder)
    [print(center) for center in centers]
