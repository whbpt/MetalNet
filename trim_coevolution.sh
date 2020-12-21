gene=`echo $1`
L=`head -n 2 $gene.msa | awk '{if(NR==2)print(length($0))}'`
val=`expr $L \* 2`
tail -n +2 $gene.pair | sort -grk 4 | head -n $val > $gene.copair 
