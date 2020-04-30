#!/bin/bash


# # Create a folder for the analysis, add the extracted fastq in a subfolder called fastq
# # and launch the analysis by typing bash pipeline_3.0.sh --fastq /full/path/to/folder/fastq --ref_dir /path/to/ref
# # Example:
# # bash /home/scripts/single_script1.5.sh  --fastq /name_of_dataset/fastq --ref-dir /data/reference_genomes



#--------------------------
#Sets up deflauts parameters.

nthreads=${nthreads:-8}
max_memory=${max_memory:-60000000000}

trim_qual_score=${trim_qual_score:-0}
trim_min_len=${trim_min_len:-60}	
trim_qual_window_size=${trim_qual_window_size:-5}

megahit_k_step=${megahit_k_step:-10}
megahit_k_min=${megahit_k_min:-21}
megahit_min_count=${megahit_min_count:-2}

concoct_length_tresh=${concoct_length_tresh:-1000}
concoct_clusters=${concoct_clusters:-400}
concoct_k_mer_length=${concoct_k_mer_length:-4}

#Named argument for multiple parameters testing
possible_params=(fastq 
ref_dir 
nthreads 
max_memory 
outdir 

trim_qual_score
trim_min_len 
trim_qual_window_size 

megahit_k_step 
megahit_k_min 
megahit_min_count 

concoct_length_tresh 
concoct_clusters 
concoct_k_mer_length)


## Overrides the default by user-specified parameters 
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
		#store the parameter name (without "--" in the variable $param 
        param="${1/--/}"
		
		#Test if $param is a modifiable parameter, other wise exits. 
		if [[ ! " ${possible_params[@]} " =~ " ${param} " ]]; then
			{ echo >&2 "Parameter: $param not recognized"; exit 1; }
		fi	 
		
		#Override the default value of $param for the user-specified one: $2
		#echo $param "$2"
        declare $param="$2"
        # echo $1 $2 // Optional to see the parameter:value result
   fi

  shift # goes to the next argument in the list (details http://tldp.org/LDP/Bash-Beginners-Guide/html/sect_09_07.html)
done


# Tests if a fastq folder was specified, otherwise exits
if [ -z $fastq ]; then
	{ echo >&2 "Must provide fastq folder"; exit 1; }
fi

# Reminds the user of the reference directory specified, if any. 
if [ ! -z $ref_dir ]; then # if the string $ref_dir (script second argument) is not empty
    echo "Reference genomes provided at:" $ref_dir
else 
	echo "No reference genomes provided"
fi

#----------------------------
# Managing inputs

# Corrects the folder path if is finishes with "/"
if [ "${fastq: -1}" == "/" ];then
    tmp=${fastq::-1}
    fastq=$tmp
fi

# If no --outdir was specified, the results will be stored at the emplacement of the fastq folder, under results/date_time/
if [ -z $outdir ]; then
	outdir=${fastq/fastq/results}
	tmp=${outdir::-1}
    outdir=$tmp
	printf -v date '%(%Y-%m-%d_%H-%M-%S)T' -1 
	echo output will be saved at $outdir/$date
	tmp=$outdir"/"$date
	outdir=$tmp
fi


#If the fastq are noted with .fq, changes it to .fastq
find $fastq -name "*.fq" -exec sh -c 'mv "$1" "${1%.fq}.fastq"' _ {} \;

# stores the name of the folder containg the fastq folder in $dataset
dataset=`echo $fastq | rev | cut -f2 -d "/" | rev `

echo The directory is $fastq
echo The dataset name is $dataset


#--------------------------
# Preparing output folder, computation time file and stores the parameter used.

if [ ! -d $outdir ] ; then
	mkdir -p $outdir 
fi

# Prepares the files to save the comuptation time of each step:
if [ ! -e $outdir/times.tab ]; then
    echo -e "Dataset"'\t'"FastQC1"'\t'"Trimmomatic"'\t'"FastQC2"'\t'"MEGAHIT(assembly)"'\t'"MetaQuast"'\t'"CONCOCT(binnig)"'\t'"CheckM(binning_qc)"'\t'"Total" > $outdir/"times.tab"
fi

echo -n $dataset > $outdir/runtimes_tmp #start row 

cat $outdir/runtimes_tmp ; echo 
total_time_start=`date +%s`

# Saves the parameters used in parameters.tab
echo "outdir: "$outdir
echo Parameters > $outdir"/parameters.tab"
echo -e "Script"'\t'"$0" > $outdir"/parameters.tab"
for i in ${possible_params[@]}; do
	echo -en "$i"'\t' >> $outdir/parameters.tab
	eval echo -e "\$$i" >> $outdir/parameters.tab
done


# # ---------------------
# # FastQC1
# Sets up output directory
raw_QC=$outdir"/raw_reads_QC"

start=`date +%s` # Stores the time at which FastQC1 started.

# if Fastqc1 wasn't already performed. 
if [[ ! -d $raw_QC ]];then
	mkdir "$raw_QC"
    
    /opt/FastQC/fastqc -o $raw_QC --noextract -f fastq --threads $nthreads  $fastq/*.fastq 
	 # -o $raw_QC: Output directory
	 # --noextract: Do not uncompress the output file after creating it.
	 # -f fastq: Bypasses the normal sequence file format detection.  
	 # --threads $nthreads: Number of threads to use
	 
	

fi
end=`date +%s`

# writes the computing time used in this step to runtimes_tmp
runtime=$((end-start))
echo -en '\t'"$runtime" >> $outdir/runtimes_tmp

# # -------------------------------------------------------
# # Trimming
# Sets up output directory
trim_out=$outdir/"trim_out"

# Sets up file for the read counts after different steps of the pipeline. 
counts=$outdir/"counts.txt"

# Sets up arrays to store the read counts
file_array=()
raw_array=()
trim_array=()

start=`date +%s` # stores the time at which the read trimming started

#requires python 2.7
source /root/miniconda3/etc/profile.d/conda.sh #to allow using "conda activate environment_name"
conda activate concoct_env

#Counts the number of reads in each fastq file before processing
for rec in `ls $fastq/*.fastq`;do file_array+=("$rec");done
for rec in `ls $fastq/*.fastq | parallel -j $nthreads -k wc -l| cut -f1 -d " " `;do raw_array+=(`echo $rec | python -c "print(int(raw_input())/4)"`);done

#Saves the results in temporary files
echo "${file_array[@]}" | tr " " $'\n' >> $outdir/tmp0.txt
echo "${raw_array[@]}" | tr " " $'\n' >> $outdir/tmp1.txt


if [[ ! -d $trim_out ]];then
	mkdir "$trim_out"
	
	for reads_R1 in $fastq/*R1.fastq;do
		filename=${reads_R1##*/}
		id=${filename%%R1*}
		reads_R2=${reads_R1/R1/R2} 
		echo launch trimming on $reads_R1 and $reads_R2
		java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $nthreads -phred33  -trimlog $trim_out/$id".log" $reads_R1 $reads_R2 $trim_out/$id"good_R1.fastq" $trim_out/$id"bad_R1.fastq" $trim_out/$id"good_R2.fastq" $trim_out/$id"bad_R2.fastq" ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:$trim_qual_window_size:$trim_qual_score MINLEN:$trim_min_len
		#  PE : paired-end mode
		# -threads $nthreads: number of threads to use
		# -phred33 : type of quality score to use 	https://en.wikipedia.org/wiki/FASTQ_format#Encoding
		# -trimlog $trim_out/$id".log" : emplacement of the summary output file
		# <inputs> <outputs>
		# ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads #adapter triming ILLUMINACLIP:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold> 
		# SLIDINGWINDOW:$trim_qual_window_size:$trim_qual_score # quality filtering: SLIDINGWINDOW<size of the window><minimum AVERAGE quality score to be kept>
		# MINLEN:$trim_min_len # size filtering: reads shorter that $trim_min_len will be discraded
         
	done 
  
fi

# Calculates the number of reads remaining after trimming
for rec in `ls $trim_out/*good_R?.fastq | parallel -j $nthreads -k wc -l | cut -f1 -d " " `;do trim_array+=(`echo $rec | python -c "print(int(raw_input())/4)"`);done

# Saves the counts arrays in the file counts.txt
echo filename > $outdir/tmp0.txt
echo "${file_array[@]}" | tr " " $'\n' >> $outdir/tmp0.txt
echo raw_count > $outdir/tmp1.txt
echo "${raw_array[@]}" | tr " " $'\n' >> $outdir/tmp1.txt
echo trimmed_count > $outdir/tmp2.txt
echo "${trim_array[@]}" | tr " " $'\n' >> $outdir/tmp2.txt
paste $outdir/tmp0.txt $outdir/tmp1.txt $outdir/tmp2.txt > $counts


end=`date +%s`

# Writes the computing time used in this step to temporary file
runtime=$((end-start))
echo -en '\t'"$runtime" >> $outdir/runtimes_tmp


# ----------------
# FastQC 2
processed_QC=$outdir"/processed_reads_QC"
start=`date +%s`

if [[ ! -d $processed_QC ]];then
	mkdir "$processed_QC"
    /opt/FastQC/fastqc -o $processed_QC --noextract -f fastq --threads $nthreads $trim_out/*good_R*.fastq

fi
end=`date +%s`

# writes the computing time used in this step to a temporary file
runtime=$((end-start))
echo -en '\t'"$runtime" >> $outdir/runtimes_tmp


# --------------------------------------
# Assembly
start=`date +%s`

assembly=$outdir"/assembly"
if [[ ! -d $assembly ]];then
	#stores path to all the trimmed reads in a variable, separated by a comma
	assembly_R1_tmp=`ls -1 $trim_out/*good_R1.fastq | tr "\n" ","`
	assembly_R1=${assembly_R1_tmp::-1}
	assembly_R2_tmp=`ls -1 $trim_out/*good_R2.fastq | tr "\n" ","`
	assembly_R2=${assembly_R2_tmp::-1}

	echo $assembly_R1 
	echo $assembly_R2 

	/opt/megahit/megahit -1 $assembly_R1 -2 $assembly_R2 -o $assembly --out-prefix megahit_out -t $nthreads -m $max_memory --k-min $megahit_k_min --k-step $megahit_k_step --min-count $megahit_min_count
fi

end=`date +%s`

# writes the computing time used in this step to a temporary file
runtime=$((end-start))
echo -en '\t'"$runtime" >> $outdir/runtimes_tmp


# --------------------------------------
# Assembly QC

start=`date +%s`

assembly_QC=$outdir"/assembly_QC"

if [[ ! -d $assembly_QC ]];then
	mkdir $assembly_QC
	
	if [ -z $ref_dir ]; then # if the string $ref_dir (script second argument) is empty, launch metaquast without reference
		python /opt/quast-5.0.2/metaquast.py  $assembly/megahit_out.contigs.fa -o $assembly_QC -t $nthreads --no-icarus --no-plots # -c contaminant? -a adapters? - map reads?
	
	else  #start metaquast with references
		ref_list_tmp=`ls -1 {$ref_dir/*.fa,$ref_dir/*.fna,$ref_dir/*.fasta}| tr "\n" ","`
		ref_list=${ref_list_tmp::-1}
		python /opt/quast-5.0.2/metaquast.py  $assembly/megahit_out.contigs.fa -o $assembly_QC -t $nthreads -r $ref_list --fragmented --no-icarus --no-plots

fi
fi
end=`date +%s`

# Writes the computing time used in this step to a temporary file
runtime=$((end-start))
echo -en '\t'"$runtime" >> $outdir/runtimes_tmp
cat $outdir/runtimes_tmp ; echo


#---------------------------------------
# Read mapping
start=`date +%s`

# Generate a folder /mapping contataing each read sample mapped bask to the original contigs
mapping=$outdir"/mapping"
if [[ ! -d $mapping ]];then
	mkdir $mapping
fi

# Indexes the assembled contigs
/opt/bwa/bwa index $assembly/megahit_out.contigs.fa

# Stores the indexing time
echo -en "index; ">$mapping/align_time.tab;echo $((`date +%s`-start)) >>$mapping/align_time.tab

# Aligns the reads to assembly
for R1 in `ls $trim_out/*good_R1.fastq`; do
    R2="${R1/R1/R2}"

    ID=${R1##*/}
    shrt_ID=${ID%%_good*}
    
    echo $R1 $R2 $shrt_ID
    /opt/bwa/bwa mem -t $nthreads $assembly/megahit_out.contigs.fa $R1 $R2 > $mapping/"$shrt_ID"".sam" 
done

# Stores the aligning time 
echo -en "map;">>$mapping/align_time.tab;echo $((`date +%s`-start)) >>$mapping/align_time.tab

# Compresses sam files to bam and sorts them
for file in `ls $mapping/*.sam`; do
  out=${file/.sam/.bam}
  echo converts $file to bam
  if [[ ! -f "$out" ]]; then
	echo $out
    samtools view -b -@ $additionnal_threads -1 $file > $out
	rm $file
  fi
  
  echo sorts $out
  memory_per_thread=$(($max_memory/$nthreads))
  additionnal_threads=$(($nthreads-1))
  samtools sort -@ $additionnal_threads -m $memory_per_thread -l 5 -o ${out/bam/sorted.bam} $out
done

# Records the time used for compression and sorting
echo -en "compress & sort; ">>$mapping/align_time.tab;echo $((`date +%s`-start)) >>$mapping/align_time.tab

# Indexes the resulting bam files
for file in `ls $mapping/*sorted.bam`; do
  while [ $(jobs -r | wc -l) -ge $nthreads ] ; do sleep 1 ; done
  samtools index $file &
done

# Waits for all to finish
while [ $(jobs -r | wc -l) -ge 1 ] ; do sleep 1 ; done

# Records the time taken to index the sorted bam files
echo -en "index2; ">>$mapping/align_time.tab;echo $((`date +%s`-start)) >>$mapping/align_time.tab

# Creates arrays to store the number of duplicates and the number of reads mapping to the assembly
assembly_array=()
duplicates_array=()
for i in $(seq 1 $N_fastq);do 
	assembly_array+=(NA)
	duplicates_array+=(NA)
done

# Stores the number of dupicates and of reads mapping to the assembly in the arrays
increment=0
for file in `ls $mapping/*sorted.bam`; do	
    R1_unmapped=`samtools view -c -@ $additionnal_threads -f 68 $file` 
	# -c : only counts
	#-f 68 : flag that codes for UNmapped R1 (https://broadinstitute.github.io/picard/explain-flags.html). 
	#Has to do it this way bacause counting mapping reads is an overstimation: some reads map to mulitple places.
    R1_trim=`echo "${trim_array[@]}" | cut -f $((increment+1)) -d " "`
	echo "Remaining after trimming in $file R1: $R1_trim"
	echo "Unmapped in $file R1: $R1_unmapped"
	echo "Mapped in $file R1: " $((R1_trim-R1_unmapped))

    assembly_array[$increment]=$((R1_trim-R1_unmapped))
	duplicates_array[$increment]=`samtools view -c -@ $additionnal_threads -f 1088 $file` #flag 1088 R1 and optical duplicate
	increment=$((increment+1))
	
    R2_unmapped=`samtools view -c -@ $additionnal_threads -f 132 $file` #-f 132 codes for unmapped R1
    R2_trim=`echo "${trim_array[@]}" | cut -f $((increment+1)) -d " "`
    assembly_array[$increment]=$((R2_trim-R2_unmapped)) 
	
	duplicates_array[$increment]=`samtools view -c -@ $additionnal_threads -f 1152 $file` #flag 1152 R1 and optical duplicate
	increment=$((increment+1))
done

# Store the content of the arrays to the counts.txt file
echo mapped > $outdir/tmp3.txt
echo "${assembly_array[@]}" | tr " " $'\n' >> $outdir/tmp3.txt
echo duplicates > $outdir/tmp4.txt
echo "${duplicates_array[@]}" | tr " " $'\n' >> $outdir/tmp4.txt 

paste $counts $outdir/tmp3.txt $outdir/tmp4.txt > $outdir/tmp
mv $outdir/tmp $counts

source /root/miniconda3/etc/profile.d/conda.sh #to allow using "conda activate environment_name"
conda activate concoct_env


# --------------------------------------
# CONCOCT
start=`date +%s`

binning=$outdir"/binning"

if [[ ! -d $binning ]];then
	mkdir "$binning"
fi
# Cut up the original contigs in 10K slices

/opt/CONCOCT-1.0.0/scripts/cut_up_fasta.py $assembly/megahit_out.contigs.fa -c 10000 -o 0 --merge_last -b $binning/contigs_10K.bed > $binning/contigs_10K.fa

# Creates the coverage table
concoct_coverage_table.py $binning/contigs_10K.bed $mapping/*.sorted.bam > $binning/coverage_table.tsv 2>$binning/coverage_table.err

# Starts concoct
concoct --composition_file $binning/contigs_10K.fa --coverage_file $binning/coverage_table.tsv -b $binning/concoct_output/ -l $concoct_length_tresh -c $concoct_clusters -k $concoct_k_mer_length

# Merges clusters results
merge_cutup_clustering.py $binning/concoct_output/clustering_gt*.csv > $binning/concoct_output/clustering_merged.csv

# Creates a folder with the bins in fasta format
mkdir $binning/concoct_output/fasta_bins
extract_fasta_bins.py $assembly/megahit_out.contigs.fa $binning/concoct_output/clustering_merged.csv --output_path $binning/concoct_output/fasta_bins

end=`date +%s`

# writes the computing time used in this step to a temporary file
runtime=$((end-start))
echo -en '\t'"$runtime" >> $outdir/runtimes_tmp
cat $outdir/runtimes_tmp ; echo

# -------------------------
# Check-M

conda activate checkm_env

b_QC=$outdir"/binning_QC"
start=`date +%s`

if [[ ! -d $b_QC ]];then
	mkdir "$b_QC"
fi

checkm lineage_wf --reduced_tree $binning/concoct_output/fasta_bins/ $b_QC/ -t $nthreads -x fa

checkm qa -o 2 $b_QC/lineage.ms $b_QC > $b_QC/extended_report.txt

end=`date +%s`

# Writes the computing time used in this step to a temporary file
runtime=$((end-start))
echo -en '\t'"$runtime" >> $outdir/runtimes_tmp
cat $outdir/runtimes_tmp ; echo

#----------
#Calculates total runtime

total_time_end=`date +%s`
total_runtime=$((total_time_end-total_time_start))

echo -e '\t'"$total_runtime" >> $outdir/runtimes_tmp

# Merges the temporary files
cat $outdir/times.tab $outdir/runtimes_tmp > $outdir/runtimes_tmp2
mv $outdir/runtimes_tmp2 $outdir/times.tab

# Removes the temporary files
rm $outdir/runtimes_tmp*
rm $outdir/tmp*
conda deactivate

