for i in `cat $1top3_model.txt`
do
raxml-ng --support --tree $1'T_'$i'.raxml.bestTree' --bs-trees $1'bs_'$i'.raxml.bootstraps' --prefix $1'S_'$i --redo --threads 4 --force
done
