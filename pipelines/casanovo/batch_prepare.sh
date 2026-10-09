#!/bin/bash

suffix=.pdf
# collapsing psms into peptide level
for file in *.tab
    do
        folder=`basename $file .tab`
        python prepare_from_casanovo.py $file $folder
    done


# counting the peptides by mapping to whole protein sequence and making files with count_per pattern
# this is the most expensive part of the code; we used ahorasick package based on indexing to make it faster
for psm in psm_*.csv
    do 
        python count.py $psm
    done

# making the visualization
for report in count_per_*.csv
    do 
        python visualize2.py -i $report -o $report$suffix
    done
