#!/bin/bash

for file in *.tab
do
folder=`basename $file .tab`
python prepare_from_casanovo.py $file $folder
done

for psm in psm_*.csv;do python count.py $psm;done
for i in `ls count_per_*.csv`;do python visualize2.py -i $i -o "$i".pdf ;done
