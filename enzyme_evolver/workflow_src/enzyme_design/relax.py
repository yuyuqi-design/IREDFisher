import subprocess as sp
relax = f'rosetta_scripts.mpi.macosclangrelease @refine.flags -nstruct 100'
sp.run(relax,shell=True)