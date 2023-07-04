#!/bin/bash

# Check alignment efficiency
############################

# Create new  file with summary
echo -n > results/data/insertSize_summary.txt

# Extract all efficiencies
selected_files=$(ls -d -1 $PWD/temp/bam_stat/*)

# Process one by one
for f in $selected_files; 
do 
	# echo 'Processing' $f;
	echo -n -e $f "\t" >> results/data/insertSize_summary.txt;
  echo -n -e $(egrep -i 'insert size average' $f | grep -o '[0-9]*\.[0-9]*') "\t" >> results/data/insertSize_summary.txt
  echo -n -e $(egrep -i 'insert size standard deviation' $f | grep -o '[0-9]*\.[0-9]*') "\n" >> results/data/insertSize_summary.txt
done
