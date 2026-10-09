#!/bin/bash
sep=/
suff=.mzml
sec=.tab
for file in *.raw
do
rm *.mztab
rm *.log
folder=`basename $file .raw`	
echo $folder$sep$folder$suff
casanovo sequence --config casanovo.yaml $folder$sep$folder$suff
for mztab in *.mztab
do
basefile=`basename $mztab .mztab`
mv $mztab $folder$sec
done



done
