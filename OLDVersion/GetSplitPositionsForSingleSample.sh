for file in $(ls *bam)
do

###Get left breakpoints for deletions

cat $file|awk '$5>=0'|awk '$NF ~ /S/'|awk 'length($NF)==5 || length($NF)==6'|\
 awk '$NF !~ /I/'|awk '$NF !~ /D/'|sed 's/M/\t/g'|\
 awk 'NF==8 && $NF ~ /S/'|awk '$(NF-1) !~ /S/'|awk '{print $3,$5}' > $file.LeftBreaks.SplitRead.txt


###Get right breakpoints for deletions

cat $file|awk '$5>=0'|awk '$NF ~ /S/'|awk 'length($NF)==5 || length($NF)==6'|\
 awk '$NF !~ /I/'|awk '$NF !~ /D/'|sed 's/S/\t/g'|\
 awk 'NF==8 && $NF ~ /M/'|awk '$(NF) !~ /S/'|awk '{print $2+1,$5}' > $file.RightBreaks.SplitRead.txt

#############All split positions

cat $file.LeftBreaks.SplitRead.txt $file.RightBreaks.SplitRead.txt|sort|uniq -c > temp; mv temp $file.Dup.SplitRead.txt

##Left split positions
cat $file.LeftBreaks.SplitRead.txt|sort|uniq -c > temp; mv temp $file.LeftBreaks.SplitRead.txt

##right split positions
cat $file.RightBreaks.SplitRead.txt|sort|uniq -c > temp; mv temp $file.RightBreaks.SplitRead.txt


done 
