awk -F, '{if(NR>1)print $1,$2}' $1 > $2
