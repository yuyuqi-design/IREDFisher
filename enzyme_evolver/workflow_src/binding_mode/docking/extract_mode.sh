#first argument is the name of protein structure name without extension
#the second argument is the name of ligand structure name without extension
#third argument is the docking_dir
#fouth argument is the job_dir
mode=`awk '{print $1}' $3scores.txt |head -1 | sed 's/log/config/;s/\.txt/_out\.pdbqt/'`
cp $3$1'.pdbqt' $4$1'_'$2'_best_mode.pdbqt'
cat $mode >> $4$1'_'$2'_best_mode.pdbqt'
#obabel $4$1'_'$2'_best_mode.pdbqt' $4$1'_'$2'_best_mode.pdb'
#sed -i '/MODEL        1/d' $4$1'_'$2'_best_mode.pdb'