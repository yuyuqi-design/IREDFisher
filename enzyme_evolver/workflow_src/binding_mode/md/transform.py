import sys
import subprocess as sp
sys.path.append('/home/yuqi/software/anaconda3/envs/ambertool16/lib/python3.6/site-packages/')
import parmed as pmd

amber = pmd.load_file('complex.prmtop', 'complex.inpcrd')
# Save a GROMACS topology and GRO file
amber.save('complex.top')
amber.save('complex.gro') 