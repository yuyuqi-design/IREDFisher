for i in `cat $2top3_model.txt`
do
raxml-ng --msa $1 --model $i'+G+F' --prefix $2'T_'$i --seed 2 --tree pars{50},rand{50} --redo --force
done
