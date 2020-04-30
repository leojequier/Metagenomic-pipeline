#!/usr/bin/env Rscript

# Reads command line aguments
args = commandArgs(trailingOnly=TRUE)

# Command:
# Rscript --vanilla Xiao_sampling.R Xiao_abundances.csv out.csv

# args=c("tmp", # abundance in, in a column format
#        "out.csv", # file name out: number of read sampled per genome/sample in csv format
#        100, # Size of the total dataset (sum of illumina lanes)
#        0.1, # Noise constraint. Between 0 and 1. 0.1, the new abundance will be + or - 10%
#        9, # Number of samples
#        300, # Read size
#        12 # Number of Illumina lanes
#        )
 
# test if there is at least 7 argument: if not, return an error
if (length(args)!=7) {
    print("argument required: abundances_in, file_out, total_size(gbp), Noise_constait,read_lenth, Nconditions")
}

# stores input in a data frame
df=read.csv(args[1], sep=';',header = T)

# Measure the number of genomes
ngenomes=nrow(df)

# Read the total size of the dataset from the agruments, convert the total size from Gbp to bp 
total_size=as.numeric(args[3])*(10**9)

# Reads the number of sample to produce from the arguments
nsamples=as.numeric(args[5])
read_size=as.numeric(args[6])

# Calculates the number of nucleotides per sample
sample_size_nuc=total_size/nsamples

# Stores the base abunance from the input file to a vector
base_abundances=df[,3]

# constrains the noise to a given percentage
noiseConstraint=as.numeric(args[4])
noiseSigma = noiseConstraint * base_abundances 

# Calculates the number of read per sample to produce
nreads=ceiling(sample_size_nuc/read_size)#number of reads/samples

# Creates a matrix to store the number of reads. Rows are samples and columns are species
m=matrix(rep(0, times=nsamples*ngenomes),ncol = ngenomes)#create a matrix to store the abundance in each samples 

#for each sample to create in the input
for(row in 1:nsamples){
	# Calculates the noise to add to the base abundance
    noise = noiseSigma * rnorm(sd = 1, n= length(noiseSigma)) #normally distributed noise mu=0 sd=1, mulitplied by the base abundance *5%
    noisy_abundance=base_abundances+noise #the noise (can be pos or neg) is added to the base abundance
    x=sample(as.factor(df$name), prob = noisy_abundance, size = nreads,replace = TRUE) #sample reads in proportion to the noisy abundance.
	# replace the matrix row with the calculated values 
    m[row,]=table(x)
}

# Cross-check to ensure the abundance levels and coverage levels obtained for C. scindes corresponds to the planned levels
csindens=which(names(table(x))=="Clostridium_scindens")
print("Abundance in each sample:")
print(m[,csindens]/nreads)

coverage=apply(X = m, MARGIN = 2, function(x){sum(x)*read_size/3500000})
print("Total coverage")
print(coverage[csindens])

abundance=apply(X = m, MARGIN = 2, function(x){sum(x)/(nreads*nsamples)})
print("Total abundance")
print(abundance[csindens])

# Converts the matrix to a data frame to facilitate exportation. Transpose it to have the species as rows.
d_export=data.frame(t(m))
# Add the species names in the first columns
d_export2=cbind(names(table(x)), d_export)
# Merge with the input data to get the reference genome ID corresponding to each species name.
d_export3=merge(df, d_export2, by.x=1, by.y =1)
# Keeps only relevant columns. 
d_export4=d_export3[,c(2,4:ncol(d_export3))]
# add the sample number as a name for each column "Sample_1; Sample_2...."
names(d_export4)[2:ncol(d_export4)]=paste(rep("Sample_", times=nsamples), 1:nsamples,sep="")
#writes the data frame to csv
write.csv(d_export4, file=args[2],row.names = F)


