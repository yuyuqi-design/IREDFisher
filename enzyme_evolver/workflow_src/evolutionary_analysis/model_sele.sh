##i=1
#for model in `cat enzyme_evolver/workflow_src/evolutionary_analysis/pro_model.txt`
#do
#raxml-ng --evaluate --msa $1 --model $model'+G+F' --tree $2'start.raxml.bestTree' --prefix $2'E_'$model --data-type AA --redo --force
##i=$((i+1))
#echo $i
#done
grep 'AIC score' $1'E'*raxml.log| awk -F / '{print $2}'|sort -n -t ':' -k3 | head -1 |awk -F : '{print $1}'| awk -F _ '{print $2}'|awk -F . '{print $1}' > $1'top3_model.txt'
    
