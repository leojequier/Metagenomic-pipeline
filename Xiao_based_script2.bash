#!/bin/bash
 
# creates a "resampled" directory and add some constrained gaussian noise to the abundance values in the input csv 
# input format: 
# name;reference;cumulated_abundance;dataset_1;...
# Achromobacter xylosoxidans;GCF_001457475;0.003401384;0.003398635; ....
# ...

# Usage: bash Xiao_based_script2.bash Xiao_abundances.csv


file=$1

#Convert the excel file to unix copatible format and replace spaces by underscores
cat $1 | dos2unix | tr " " "_" > spaceless.csv

# Measure the number of columns of the input csv
ncol=`head -1 spaceless.csv | tr ";" " " | wc -w`
echo $ncol

# If the output directory does not exist, creates it.
if [ ! -d  resampled ] ; then 
	mkdir resampled
fi

# For each abundanc level in the input csv
for abundance_lvl in `seq 4 $ncol`;  do
	# Stores the output file name in a variable
	out=`head -1 spaceless.csv | cut -f $abundance_lvl -d";"`
	
	# Creates a temporary file respecting the input requirement of the r script
	cat spaceless.csv | cut -f 1,2,$abundance_lvl -d";" > rscript_input.tmp

	# Launch the R script on the csv file column
	Rscript --vanilla Xiao_sampling.R rscript_input.tmp "resampled/""$out""_resampled.csv" 100 0.1 9 300 12
	# arguments: 
		# [1] (tmp) the input file for the r script. Has 3 columns: Reference genome ID; Base abundance; Target Abundance
		# [2] ("resampled3/""$out""_resampled.csv") output file name
		# [3] (100) size of the dataset, in Gbp
		# [4] (0.1) Constraint to the gaussian noise, ]0;1]. 0.1 => the new abundance will be between [base_abundance-10%;base_abundance+10%]
		# [5] (9) number of samples to creates
		# [6] (300) Read length
		# [7] (12) number of illumina lanes used (useless in final design, leave at 12)
		
	
done
