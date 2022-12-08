#run as parse mode to generate binary alignmnet file (*.raxml.rba)
#raxml-ng-mpi --parse --msa $1 --prefix parse --model LG  --seed 2
#using RBA alignment to generate the starting tree (*.raxml.bestTree)
raxml-ng --msa $1 --prefix start --model LG  --seed 2
#check the reasonability of the starting tree and evaulate the protein substitution model and output the top3_model.txt
sh ../../workflow_src/evolutionary_analysis/model_sele.sh $1
#compute the tree and calculate if the tree is the same using rfdistance
sh ../../workflow_src/evolutionary_analysis/top3model_tree.sh $1
#bootstrap for 3 trees
sh ../../workflow_src/evolutionary_analysis/top3model_tree_bs.sh $1
#supported values
sh ../../workflow_src/evolutionary_analysis/top3model_tree_support.sh $1
sh ../../workflow_src/evolutionary_analysis/clean.sh