#!/usr/bin/env python3
#Author: Dr. Yuqi YU. yuqi.yu@manchester.ac.uk.
#This script is used to make constraint for PROTEIN,Heavy,MC,CA and generate posre and modify the topol or itp file for including the restraint file
#warning: only apply to protein with only one continuous chain
#usage: ./make_restraint.py proteinStructure.pdb(recommend pdb2gmx_process.pdb) topol.top(or itp file)

import subprocess as sp   
import sys

def gmx_restraint(input='complex.gro',folder=''):
    
    #command
     
    restr_Protein_comm = 'echo 1|gmx genrestr -f ' + str(folder+input) + ' -o ' + folder+'posre_PROTEIN.itp'
    restr_Heavy_comm = 'echo 2|gmx genrestr -f ' + str(folder+input) + ' -o ' +folder + 'posre_Heavy.itp'
    restr_MC_comm = 'echo 5|gmx genrestr -f ' + str(folder+input) + ' -o ' + folder+'posre_MC.itp'
    restr_CA_comm = 'echo 3|gmx genrestr -f ' + str(folder+input) + ' -o ' + folder+'posre_CA.itp'
    
    #call the commmand genrestr to generate restraint file
    
    print('starting to generate protein restraint...\n') 
    sp.run(restr_Protein_comm, shell = True)
    print('starting to generate Heavy atoms restraint...\n')
    sp.run(restr_Heavy_comm, shell = True)
    print('starting to generate MC atoms restraint...\n')
    sp.run(restr_MC_comm, shell = True)
    print('starting to generate CA atoms restraint...\n')
    sp.run(restr_CA_comm, shell = True)
    sp.run('rm \#*', shell = True)

if __name__ == '__main__':
        input=sys.argv[1]
        gmx_restraint(input)

# ###############################################################################################
# ###############################################################################################
#
# input = sys.argv[1]
# gmx_restraint(input)
#
# file = sys.argv[2]
# restraint_to_itp(file)
