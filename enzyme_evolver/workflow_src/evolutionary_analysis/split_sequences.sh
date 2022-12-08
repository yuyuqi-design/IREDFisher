grep -F  '>' $1 |awk -F \> '{print $2}'|while read i; do grep $i -A 1 $1 > $i'.fasta'; done
