#Usage: sh N_rename.sh input_file N_index output_file
awk -v N_index=$2 '{if(NR<N_index)print}' $1 > $3
awk -v N_index=$2 '{if(NR==N_index) print}' $1 |sed -E 's/(.{13}).{3}/\1N01/' >>$3
awk -v N_index=$2 '{if(NR>N_index)print}' $1 >>$3
