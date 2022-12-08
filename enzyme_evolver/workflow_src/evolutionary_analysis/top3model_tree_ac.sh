for i in `cat $2top3_model.txt`
do
raxml-ng --msa $1 --model $i'+G+F'  --tree $2'T_'$i'.raxml.bestTree' --ancestral --prefix $2'anc_'$i --data-type AA --redo --force --threads 4
done
