for i in `cat $2top3_model.txt`
do
raxml-ng --bootstrap --bs-trees 100 --msa $1 --model $i'+G+F' --prefix $2'bs_'$i --seed 2 --tree pars{20},rand{20} --redo --threads 12 --force
done
