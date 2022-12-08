gmx_mpi grompp -c complex.gro -f ./min_1.mdp -o min1.tpr -p topol_restraint.top
mpirun -np 96 gmx_mpi mdrun -v -deffnm min1
gmx_mpi grompp -c min1.gro -f ./min_2.mdp -o min2.tpr -p topol_restraint.top
mpirun -np 96 gmx_mpi mdrun -v -deffnm min2
gmx_mpi grompp -c min2.gro -f ./min_3.mdp -o min3.tpr -p topol_restraint.top
mpirun -np 96 gmx_mpi mdrun -v -deffnm min3
gmx_mpi grompp -c min3.gro -f ./min_4.mdp -o min4.tpr -p topol_restraint.top
mpirun -np 96 gmx_mpi mdrun -v -deffnm min4
gmx_mpi grompp -c min4.gro -f ./nvt.mdp -o nvt.tpr -p topol_restraint.top -maxwarn 2
mpirun -np 96 gmx_mpi mdrun -v -deffnm nvt
gmx_mpi grompp -c nvt.gro -f ./npt.mdp -o npt.tpr -p topol_restraint.top -maxwarn 2
mpirun -np 96 gmx_mpi mdrun -v -deffnm npt
gmx_mpi grompp -c npt.gro -f ./npt2.mdp -o npt2.tpr -p topol_restraint.top -maxwarn 2
mpirun -np 96 gmx_mpi mdrun -v -deffnm npt2
gmx_mpi grompp -c npt2.gro -f ./npt3.mdp -o npt3.tpr -p topol_restraint.top -maxwarn 2
mpirun -np 96 gmx_mpi mdrun -v -deffnm npt3
gmx_mpi grompp -c npt3.gro -f ./npt4.mdp -o npt4.tpr -p topol_restraint.top -maxwarn 2
mpirun -np 96 gmx_mpi mdrun -v -deffnm npt4
gmx_mpi grompp -c npt4.gro -f ./npt5.mdp -o npt5.tpr -p topol_restraint.top -maxwarn 2
mpirun -np 96 gmx_mpi mdrun -v -deffnm npt5
gmx_mpi grompp -c  npt5.gro -f ./pr1.mdp -o  pr1.tpr  -p topol_restraint.top -maxwarn 2
mpirun -np 96 gmx_mpi mdrun -deffnm pr1 

#gmx mdrun -deffnm min1 -ntomp 16
