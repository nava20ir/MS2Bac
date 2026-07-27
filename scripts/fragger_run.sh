#!/bin/bash
# usage: script.sh raw_file_name
raw_file_name=`basename $1 .raw`
suffix=.tsv
suffix_param=.params
for fasta in *.fasta
do
	echo running $fasta
	python make_fragger_parameter_file.py $fasta $fasta$suffix_param
	docker run -it -v "$PWD/:/fragger_projects"  fragpipebbm java -Xmx32G  -jar /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/MSFragger-4.4.1/MSFragger-4.4.1/MSFragger-4.4.1.jar /fragger_projects/$fasta$suffix_param  /fragger_projects/$1
	mv $raw_file_name$suffix $raw_file_name$fasta$suffix
done