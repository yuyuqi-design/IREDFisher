#using alignment to generate the starting tree (*.raxml.bestTree)
raxml-ng --msa $1 --prefix $2start --model LG  --seed 2 --tree pars{10},rand{10} --redo --threads 4 --force
#check the reasonability of the starting tree and evaulate the protein substitution model and output the top3_model.txt
sh enzyme_evolver/workflow_src/evolutionary_analysis/model_sele.sh $1 $2
#compute the tree and calculate if the tree is the same using rfdistance
sh enzyme_evolver/workflow_src/evolutionary_analysis/top3model_tree.sh $1 $2
#bootstrap for 3 trees
sh enzyme_evolver/workflow_src/evolutionary_analysis/top3model_tree_bs.sh $1 $2
#supported values
sh enzyme_evolver/workflow_src/evolutionary_analysis/top3model_tree_support.sh $2
#ancestral sequence calculation
sh enzyme_evolver/workflow_src/evolutionary_analysis/top3model_tree_ac.sh $1 $2
#clean the directory
sh enzyme_evolver/workflow_src/evolutionary_analysis/clean.sh $2

