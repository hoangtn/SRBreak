
#!/bin/bash

line=Chr21.simulateDeletionsFrom1000Genomes.fa

mkdir PairedEnd
mkdir SingleEnd

temp=${line}

j=48162246

for i in 1 2 3 4 5
do

##Read length = 100; divide only 100
reads=$(expr $(expr $j \* $i) \/ 200)

dwgsim -C $i -1 100 -2 100 -c 0 -H $line PairedEnd/$temp.$i.x

##Read length = 100; divide only 100
reads=$(expr $(expr $j \* $i) \/ 100)

dwgsim -C $i -1 100 -2 0 -c 0 -H $line SingleEnd/$temp.$i.x


done



