## Abstract
Bile acids (BA) transforming bacteria, or more specifically, bacteria performing 7-dehydroxylation (7-DH), are suspected to play a crucial part in the endocrine role of BA and the interaction between the gut microbiota and the host. 7-DH bacteria’s low abundance levels in the gut microbiota, however, render them difficult to study. During this project, a bioinformatic pipeline for the analysis of whole metagenomic shotgun sequencing datasets was developed with the aim of maximizing sensitivity towards low abundance species, thus enabling the potential discovery of new 7-DH bacteria. This pipeline was consolidated into a singular script, which enables parameters tested in this study to be easily modified by command line arguments, facilitating the usage for future users. Further, this pipeline has been containerized with docker so that specific versions of tools and dependencies are pre-installed, making the results reproducible and the pipeline compatible across computers. Finally, realistic artificial datasets simulating the mouse gut microbiota were created to optimise and validate the pipeline. With our validated parameters, the pipeline produced a Metagenome-Associated Genome (MAG) with 80% completeness and <5% contamination for C. scindens, a model for 7-DH bacteria, at an abundance of 0.05%. These results indicate that for most of the abundance range at which 7-DH bacteria are expected, the pipeline presented here will likely be able to produce MAGs with a sufficient quality to better grasp the diversity of 7-DH bacteria and discover new genes related to BA metabolism.

For more details, see MS thesis PDF document.

## Installation: 
1)Download Docker image:
...
2)Create container
..

## Usage : 

/opt/home/scritps/pipeline_3.0.sh --fastq <fastq_folder> [ --ref_dir --nthreads --max_memory --outdir --trim_qual_score --trim_min_len --trim_qual_window_size --megahit_k_step --megahit_k_min --megahit_k_max --megahit_min_count --concoct_length_tresh --concoct_clusters --concoct_k_mer_length]


