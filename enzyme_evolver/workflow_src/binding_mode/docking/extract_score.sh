for i in `ls -1 $1log*.txt`;do echo -n $i'  ';grep -A1 '\----' $i|tail -1|awk '{sum+=$2}END{print sum/NR}'; done|sort -n -k2> $1'scores.txt'
#