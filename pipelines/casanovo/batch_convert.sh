#!/bin/bash
for file in *.raw
do
folder=`basename $file .raw`
python raw2mzml.py -i $file -o $folder
done
