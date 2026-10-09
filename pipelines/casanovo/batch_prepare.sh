#!/bin/bash


# collapsing psms into peptide level
for file in *.tab
    do
        folder=`basename $file .tab`
        python prepare_from_casanovo.py $file $folder
    done


# counting the peptides by mapping to whole protein sequence
for psm in psm_*.csv
    do 
        python count.py $psm
    done

# making the visualization
for report in `ls count_per_*.csv`
    do 
        python visualize2.py -i $report -o "$i".pdf 
    done
