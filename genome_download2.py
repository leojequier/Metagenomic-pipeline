from Bio import Entrez
import urllib.request	

Entrez.email = "leonard.jequier@unil.ch"

ref_dict={}

with open("ref_id_out.tab", "w") as data_out: #opens a file to write the reference genome found
	with open("taxid_report_clean.tab", "r") as data_in: # read the file with the taxononmy information for each species. With the form: Species_name\tspecies_taxid genus_taxid family_taxid order_taxid ...
		
		for line in data_in: #for each species

			sep=line.split("\t")
			species_name=sep[0]
			
			lowest_leaf=sep[1].split(" ")[0] # stores the the lowest taxonomic leaf (species or genus) in a variable
			
			# Searches for a reference-grade assembly on the NCBI assembly database for that species.
			handle = Entrez.esearch(db="assembly", term='txid{0}[Organism] AND "reference genome"[filter]'.format(lowest_leaf),sort="Contig N50 (down)")
			record = Entrez.read(handle)
			if len(record['IdList']) > 0: # If at least one reference-grade assembly is available
				print("Has a reference-grade genome")
				rec_id=record['IdList'][0] # Selects the first one, as it is sorted by decreasing N50 values
				
				# Extracts the FTP path to the refseq database from the text summary file.
				summary_h=Entrez.esummary(db="assembly",id=rec_id, retmode="txt")
				summary=Entrez.read(summary_h)
				ref_seq_path=summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq'] 
				id=ref_seq_path.split("/")[-1]
				
				# Downlads the reference genome selected under the folder ./genomes/
				urllib.request.urlretrieve("{0}/{1}_genomic.fna.gz".format(ref_seq_path,id), "genomes/{0}_genomic.fna.gz".format(id))
				
			#else if no reference-grade assembly can be found
			else:
				# Searches for a representative-grade assembly on the NCBI assembly database for that species.
				handle = Entrez.esearch(db="assembly", term='txid{0}[Organism] AND "representative genome"[filter]'.format(lowest_leaf),sort="Contig N50 (down)")
				record = Entrez.read(handle)
				if len(record['IdList']) > 0: #If at least one representative-grade assembly is found
					print("Has a representative-grade genome")
					rec_id=record['IdList'][0] # Selects the first one, as it is sorted by decreasing N50 values
					
					# Extracts the FTP path to the refseq database from the text summary file.
					summary_h=Entrez.esummary(db="assembly",id=rec_id, retmode="txt")
					summary=Entrez.read(summary_h)
					ref_seq_path=summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
					id=ref_seq_path.split("/")[-1]
					
					# Downlads the reference genome selected under the folder ./genomes/
					urllib.request.urlretrieve("{0}/{1}_genomic.fna.gz".format(ref_seq_path,id), "genomes/{0}_genomic.fna.gz".format(id))
				
				#else if no representative-grade assembly can be found
				else: 
					# Searches for a complete-grade assembly on the NCBI assembly database for that species.
					handle = Entrez.esearch(db="assembly", term='txid{0}[Organism] NOT partial[filter]'.format(lowest_leaf),sort="Contig N50 (down)")
					record = Entrez.read(handle)
					if len(record['IdList']) > 0:  #If at least one complete-grade assembly is found
						print("Has complete genome")
						rec_id=record['IdList'][0] # Selects the first one, as it is sorted by decreasing N50 values
						
						# Extracts the FTP path to the refseq database from the text summary file.
						summary_h=Entrez.esummary(db="assembly",id=rec_id, retmode="txt")
						summary=Entrez.read(summary_h)
						ref_seq_path=summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
						id=ref_seq_path.split("/")[-1]
						
						
						if ref_seq_path != "": # If a refseq path is available
							print("Refseq path found")
							# Downlads the reference genome selected under the folder ./genomes/
							urllib.request.urlretrieve("{0}/{1}_genomic.fna.gz".format(ref_seq_path,id), "genomes/{0}_genomic.fna.gz".format(id))
						else:
							# print an error message and indicate in the output that the geneome will hat to be manually selected
							print("{0} will have to be manually curated".format(species_name))
							id="manual curation"
							
					else:
						# print an error message and indicate in the output that the geneome will hat to be manually selected
						print("{0} will have to be manually curated".format(species_name))
						id="manual curation"
			
			 # parse the Refseq path to get the file name
			new_id="".join(id.split(".")[:-1])
			# writes the the output file the species name and the file name corresponding reference genome 
			data_out.write("{0}\t{1}\n".format(species_name, new_id))
			
			
			
