for code in `cat $1codes.txt`;
do
#cd $i;
best_file=`tail -10 $1$code/log |grep 'pdb '| sort -n -k2| awk '{print $1}' |head -1`
echo $best_file
cp $1$code/$best_file $1
#cd ../
done
