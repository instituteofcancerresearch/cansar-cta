import os
import bisect
import argparse
import csv

# This script is to create a temporary mutation file which can be parsed easier compared to pure vcf.
# It takes vcf as input and produces meta file with different mutations separated with commas.
# letmeknow=True lines are introduced for debugging purposes

# Require three inputs (filenames), recieved from command line args
#
# id_mapping_file
#   A file generate by snpEff (gene.txt by default)
#
# annotation_file
#   A VCF file annotated with snpEff and snpSift
#
# maf_file
#   a four-column, space separated file containing:
#	 the mutation index number, seq-project, patient-id and mutation location.


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--id_mapping_file')
	parser.add_argument('-a', '--annotation_file')
	parser.add_argument('-m', '--maf_file')
	args = parser.parse_args()
	id_mapping_file = args.id_mapping_file
	annotation_file = args.annotation_file
	maf_file = args.maf_file


class Mutation:
	def __init__(self):
		self.mutation_chr_pos = 0
		self.transcripts=[]
		self.mutation_pred_rankscore_of_canonical=0
		self.mutation_phastcons_of_canonical=0
		self.aa_change_singleletter = ""
		self.mutation_depth = 0

class Gene:
	def __init__(self):
		self.genename=""
		self.ensg_id=""
		self.canonical_transcript_id=""
		self.mutposlist=[]
		self.mutations=[]
		self.mutation_patient=[]
		self.mutation_study=[]

class Transcript:
	def __init__(self):
		self.enstid=""
		self.mutation_type=""
		self.mutation_impact=""
		self.mutation_protein_pos=0
		self.aa_change = ""


# Create a dict to store  ENST IDs as keys associated
# with ENGS IDs as values
id_dict = {}
with open (id_mapping_file, "r") as f:
	reader = csv.reader(filter(lambda row: row[0]!='#', f), delimiter = "\t")
	for r in reader:
		# r[1] should be ENSG ID and r[2] should be ENST ID
		# This initially confused me as I thought we would be mapping
		# multiple transcripts as keys to individual genes as values
		# Not the case, keys are genes and canonical transcripts are values
		if (not r[1].startswith("ENSG") or not r[2].startswith("ENST")):
			raise Exception('When reading the snpEff_genes.txt file, expected ENSG and ENST strings at the start of the second and thrird columns ')
		id_dict[r[1]] = r[2]

# This should be looked at later - currently just a list (?) of
# strings corresponding to lines in the file. These are later
# split in various ways and used. Wonder if three dicts indexed
# by the mutation index would be clearer?
file_maf=open(maf_file,"r")
maf_lines=file_maf.readlines()

# Used to keep track of genes that have been seen before.
ensemblgnamelist=[]

# Holds Gene objects created below
genelist=[]

#Do not consider these mutations as these are either non-coding or not well defined
unwanted_mutation_types=["TF_binding_site_variant","protein_protein_contact","TFBS_ablation","intron_variant","intragenic_variant","structural_interaction_variant","non_coding_transcript_exon_variant"]

#READING THE MUTATION VCF FILE and LINKING TO MAF FILE INFORMATION
with open(annotation_file,"r") as f:
	for line in f:
		if not line.startswith("#"): # skip headers
			line=line.rstrip("\n")
			
			#initialize some of the variables that will change with every line
			aachangesingleletter = ""
			mutation_rankscore=0
			mutation_phastcons=0
			
			# New Mutation object
			newmutation=Mutation()
			
			#line_parser=line.split(";") # Is this used anywhere?
			
			vcf_fields = line.split("\t")
			info_fields = vcf_fields[7].split(";")
			
			# [JAMES] debug print
			#print("working on: %s" % (vcf_fields[2]))
			
			# Check each field in INFO (sep by ;) for things we need
			for info_field in info_fields: # Parse annotation field of the line
				# for loop above was 'for string in line.split(";"):'
				if info_field.startswith("dbNSFP_MetaSVM_rankscore"): #Rankscore=Mutation Impact Prediction
					temp_score=info_field.split("=")[1].split(",")
					for score in temp_score:
						# [JAMES] Need to check what the comma-sep dbNSFP_MetaSVM_rankscore values
						# mean. Why are there more than one? In the code below, we accept the first
						# one that is not '.'
						if score!=".":
							mutation_rankscore=score#string.split("=")[1]
							break
						else:
							# Let me know will be True if *any* of the dbNSFP_MetaSVM_rankscore values are '.'
							letmeknow=True
				elif info_field.startswith("dbNSFP_phastCons100way_vertebrate"): # Get PhastCons Conservation Score
					temp_phastcons=info_field.split("=")[1].split(",")
					for phastcons in temp_phastcons:
						if phastcons!=".":
							mutation_phastcons=phastcons#string.split("=")[1]
							break
						else:
							letmeknow=True
				elif info_field.startswith("dbNSFP_aaref"): #Reference Aminoacid
					# [JAMES] It is not clear to me why the strings are concatenated
					# in the line below and the similar lines in the remaining elif blocks
					aachangesingleletter=aachangesingleletter+info_field.split("=")[1]
				# These two were not requested from snpEff so do not exist
				#elif info_field.startswith("dbNSFP_aapos"): # Aminoacid Position of Mutation
				#	aachangesingleletter=aachangesingleletter+info_field.split("=")[1]
				#elif info_field.startswith("dbNSFP_aaalt"): #Aminoacid result of Mutation
				#	aachangesingleletter=aachangesingleletter+info_field.split("=")[1]
				elif info_field.startswith("ANN"): #annotations for a set of comma-sep transcripts?
					transcripts = info_field.split("=")[1].split(",")
			
			# [JAMES] debug print
			#print("  After parsing INFO:\n  rank score:%s\n  phastcons: %s" % (mutation_rankscore, mutation_phastcons))
			
			# newmutation is an object of class Mutation
			newmutation.aa_change_singleletter=aachangesingleletter # Not requested from snpEff - should we add it or get if from the change for each transcript?
			newmutation.mutation_pred_rankscore_of_canonical = mutation_rankscore
			newmutation.mutation_phastcons_of_canonical = mutation_phastcons
			
			#Parse Different Transcript information for same Mutation
			#transcripts=line.split(";")[0].split(",")
			# commented out assignment of transcripts here - did it when parsing info_field
			#
			# I think the line above may rely on the first ';' appearing
			# before the dbNSFP annotations and the first ',' appearing 
			# between annotations. However, I think the ',' values separate
			# different transcript's variants
			
			# [JAMES] debug print
			#print("length of transcripts: %d" % (len(transcripts)))
			
			# split the id field in the vcf into the index number and variant allele fraction (known here as depth)
			mutation_index = int(vcf_fields[2].split("_")[0])
			mutation_depth = float(vcf_fields[2].split("_")[1])
			
			newmutation.mutation_depth=mutation_depth
			
			# Grab the line number from maf corresponding to the
			# VCF row mutation index
			# Need to add + 1 to mutation index as mutations are
			# numbered from 1
			maf_fields = maf_lines[mutation_index - 1].split(" ")
			
			# Check that the row of data we are using from the maf file matches
			# the id for the mutation from the VCF id field.
			if mutation_index != int(maf_fields[0]):
				raise Exception('There is probably an error in the .maf row ordering - mutation_index should match the first item in each row by row number')
			
			# assign study, patient and mutation key (chr_pos)
			study_id = maf_fields[1] # Get study information from Maf file
			patient_id = maf_fields[2] # TCGA patient ID
			mutation_chr_pos = maf_fields[3].rstrip() # Chromosome position of mutation
			
			
			for transcript in transcripts:
				# split this transcript's annotations using '|'
				transcript_fields = transcript.split("|")
				# Using:
				# transcript_fields[1] mutation consequence
				# transcript_fields[2] impact
				# transcript_fields[3] Gene name
				# transcript_fields[4] ENSG - Comment suggested this should be ENST
				# http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf suggestes [4] should be 'transcript' (feature type)
				# transcript_fields[6] should be ENST ID
				# transcript_fields[10] aa change
				# transcript_fields[12] position (CDS position according to docs) (element [0] after splitting on '/') 
				# transcripts_fields[13] should be the protein position
				
				study_id = maf_lines[mutation_index].split(" ")[1] #Get study information from Maf file
				patient_id = maf_lines[mutation_index].split(" ")[2].rstrip() #TCGA patient ID
				mutation_chr_pos = maf_lines[mutation_index].split(" ")[3].rstrip() #Chromosome position of mutation
				
				# [JAMES] debug print
				#print("    For the current transcript:\n    study: %s\n    patient: %s\n    chr_pos:%s" % (study_id,patient_id,mutation_chr_pos))
				
				# potential problem with next line - consequences can be concatenated with '&'
				# these won't be in the list - would be better to search for a pattern of consequences
				# that we do want to include....
				if (transcript_fields[1] not in unwanted_mutation_types) and (transcript_fields[2] != "MODIFIER"):
					
					genename = transcript_fields[3]
					ensembl_id = transcript_fields[4]
					
					newtranscript = Transcript()
					newtranscript.enstid = transcript_fields[6].split('.')[0]
					newtranscript.mutation_type = transcript_fields[1]
					newtranscript.mutation_impact = transcript_fields[2]
					newtranscript.mutation_protein_pos = transcript_fields[13].split("/")[0]
					newtranscript.aa_change = transcript_fields[10]
					
					# New transcript will get appeneded to newmutation if not seen before
					
					# [JAMES] debug print
					#print("      Using this transcript...")
					#print("      transcript ID: %s\n      consequence: %s\n      impact: %s\n      prot pos: %s\n      AA: %s" % (transcript_fields[6],transcript_fields[1],transcript_fields[2],transcript_fields[13].split("/")[0], transcript_fields[10]))
					
					#add mutation information of this transcript to the mutation class
					newmutation.transcripts.append(newtranscript)
					newmutation.mutation_chr_pos = mutation_chr_pos
					if ensembl_id not in ensemblgnamelist: # if mutation in this Gene is not added before create new Gene instance
						# [JAMES] debug print
						#print("Adding "+ensembl_id+" ("+genename+") to the genelist")
						newgene = Gene()
						newgene.genename = genename
						newgene.ensg_id = ensembl_id
						newgene.mutations.append(newmutation)
						newgene.mutation_patient.append(patient_id)
						newgene.mutation_study.append(study_id)
						newgene.mutposlist.append(mutation_chr_pos)
						# ??? is the ENSG ID (in ensembl_id) supposed to be a key in id_dict????
						# That would mean only one transcript can be stored as a key...
						# Ahhhh - id_dict stores canonoical transcript for each gene...
						if ensembl_id in id_dict.keys():
							newgene.canonical_transcript_id=id_dict[ensembl_id] #find canonical transcipt through the dictionary file that was defined in the beginning
						else:
							# This should never happen and should probably raise an exception.
							print(ensembl_id) #if ensembl ID is not found then print those IDs.
						genelist.append(newgene) #add the newgene to List of Genes
						ensemblgnamelist.append(ensembl_id) # add Gene name to List of Genenames for indexing purposes
					else: 
						# there is already mutation in this gene
						# find index of the gene so that new mutation can be added to Gene
						genename_index = ensemblgnamelist.index(ensembl_id) 
						# Check will be True if the chr_pos has not already been defined for this gene
						# or
						update_genelist = False
						if mutation_chr_pos not in genelist[genename_index].mutposlist: # if mutation position is already added before then Check=True
							update_genelist=True
						else: # if the new mutation location is not already added before
							#posindex=genelist[genename_index].mutposlist.index(mutation_chr_pos)
							# [JAMES] Not sure what the next line is supposed to do...
							#   Breaking posindex down:
							#     length of the mutposlist  minus  the index position of mutation_chr_pos counting from the back of the list  minus 1
							posindex = len(genelist[genename_index].mutposlist) - genelist[genename_index].mutposlist[::-1].index(mutation_chr_pos) - 1
							if patient_id != genelist[genename_index].mutation_patient[posindex]: 
								#if mutation position is already there but if mutation belongs to another patient then make Check=True
								update_genelist=True
							else:
								# unnecessary because check is False by default?
								update_genelist=False
						if update_genelist==True: #if this mutation is not added before per checks in previous stage
							genelist[genename_index].mutations.append(newmutation)
							genelist[genename_index].mutation_patient.append(patient_id)
							genelist[genename_index].mutation_study.append(study_id)
							genelist[genename_index].mutposlist.append(mutation_chr_pos)



#PARSING FINISHED. BELOW IS FOR CREATING OUTPUT FILES

output_file=open("annotation_summary_output_20191026.txt","w") #FILE NAME FOR OUTPUT

#PRINT HEADER FIRST - Need to check if "Mutations" is ever output as a column...
output_file.write("EnsembGenomeID"+"\t"+"Gene Name"+"\t"+"EnsemblTranscriptID"+"\t"+"\t"+"Patients"+"\t"+"Studies"+"\t"+"Pred_Rankscore"+"\t"+"PhastCons"+"\t"+"Chr Pos"+"\t"+"Mutation Type"+"\t"+"Protein Pos"+"\t"+"Mutation Impact"+"\t"+"Mutation Depth")
output_file.write("\n")

#PRINT THE METAFILE FOR EACH GENE
for gene in genelist:
	
	# [JAMES] debug print
	print("Doing a new gene object")
	
	# mutnottoadd seems to hold a list of index positions that
	# that we need to skip if the ENST ID is not the canonical
	mutnottoadd=[]
	mutindex = 0
	# Goal below is to have a list of index positions to ignore for non-canonical transcripts
	for mut in gene.mutations:
		addtooutput = False
		for enst in mut.transcripts:
			# [JAMES] debug print
			#print("mutindex: %d" % (mutindex))
			# So most things do not get added? ... Think there is a problem with recognising canonical transcripts?
			# Yuh-huh - Need to have consistency with the .version number:
			#    This Tscrpt is ENST00000369915.7, Canon Tscrpt is ENST00000369915
			
			# [JAMES] debug print
			print("This Tscrpt is %s, Canon Tscrpt is %s" % (enst.enstid, gene.canonical_transcript_id))
			if enst.enstid == gene.canonical_transcript_id:
				# [JAMES] debug print
				print("transcript is canonical")
				addtooutput=True
		if addtooutput==False:
			mutnottoadd.append(mutindex)
		mutindex = mutindex + 1
	
	output_file.write(gene.ensg_id)
	output_file.write("\t")
	output_file.write(gene.genename)
	output_file.write("\t")
	output_file.write(gene.canonical_transcript_id)
	output_file.write("\t")
	
	# b is an index similar to mutindex and checked against mutnottoadd
	b = 0
	for patient in gene.mutation_patient:
		if b not in mutnottoadd:
			output_file.write(patient)
			# This if checks if we should add a comma because
			# there are more csv elements to add
			if b < len(gene.mutation_patient)-1:
				output_file.write(",")
		b = b + 1
	output_file.write("\t")

	b = 0
	for study in gene.mutation_study:
		if b not in mutnottoadd:
			output_file.write(study)
			if b < len(gene.mutation_study)-1:
				output_file.write(",")
		b = b + 1
	output_file.write("\t")

	b=0
	for mut in gene.mutations:
		if b not in mutnottoadd:
			output_file.write(str(mut.mutation_pred_rankscore_of_canonical))
			if b<len(gene.mutations)-1:
				output_file.write(",")
		b=b+1
	output_file.write("\t")

	b = 0
	for mut in gene.mutations:
		if b not in mutnottoadd:
			output_file.write(str(mut.mutation_phastcons_of_canonical))
			if b<len(gene.mutations)-1:
				output_file.write(",")
		b=b+1
	output_file.write("\t")

	b = 0
	for mut in gene.mutations:
		if b not in mutnottoadd:
			output_file.write(str(mut.mutation_chr_pos))
			if b < len(gene.mutations)-1:
				output_file.write(",")
		b = b + 1
	output_file.write("\t")

	b = 0
	for mut in gene.mutations:
		for enst in mut.transcripts:
			if enst.enstid == gene.canonical_transcript_id:
				output_file.write(str(enst.mutation_type))
				break

		if b < len(gene.mutations)-1:
			output_file.write(",")
		b = b + 1
	output_file.write("\t")

	b = 0
	for mut in gene.mutations:
		for enst in mut.transcripts:
			if enst.enstid==gene.canonical_transcript_id:
				output_file.write(str(enst.mutation_protein_pos))
				break

		if b < len(gene.mutations)-1:
			output_file.write(",")
		b = b + 1

	output_file.write("\t")

	b = 0
	for mut in gene.mutations:
		for enst in mut.transcripts:
			if enst.enstid == gene.canonical_transcript_id:
				output_file.write(str(enst.mutation_impact))
				break

		if b < len(gene.mutations)-1:
			output_file.write(",")
		b = b + 1

	output_file.write("\t")

	b = 0
	for mut in gene.mutations:
		output_file.write(str(mut.mutation_depth))
		if b < len(gene.mutations)-1:
			output_file.write(",")
		b = b + 1


	output_file.write("\n")


output_file.close()

#OUTPUT FILE FOR LOLLIPOP PLOT, WHICH HAS ALL MUTATIONS FOR ALL TRANSCRIPTS
output_file=open("file_for_Mutation_LollipopPlot_output_20191026.txt","w")

#PRINT HEADER FILE FIRST
output_file.write("EnsemblGenomeID"+"\t"+"GeneName"+"\t"+"Canonical Transcript ID"+"\t"+"StudyName"+"\t"+"PatientName"+"\t"+"Ensembl Transcript ID"+"\t"+"MutationPos"+"\t"+"MutationType"+"\t"+"Severity"+"\t"+"Aminoacid Change"+"\t"+"Prediction rankscore-MetaSVM"+"\t"+"PhastCons Conservation Score")
output_file.write("\n")

#PRINT THE PER GENE INFORMATION
for gene in genelist: # OUTPUT FOR LOLLIPOP PLOTS
	b = 0
	for mutation in gene.mutations:
		for transcript in mutation.transcripts:
			output_file.write(gene.ensg_id)
			output_file.write("\t")
			output_file.write(gene.genename)
			output_file.write("\t")
			output_file.write(gene.canonical_transcript_id)
			output_file.write("\t")
			output_file.write(gene.mutation_study[b])
			output_file.write("\t")
			output_file.write(gene.mutation_patient[b])
			output_file.write("\t")
			output_file.write(transcript.enstid)
			output_file.write("\t")
			output_file.write(str(transcript.mutation_protein_pos))
			output_file.write("\t")
			output_file.write(transcript.mutation_type)
			output_file.write("\t")
			output_file.write(transcript.mutation_impact)
			output_file.write("\t")
			output_file.write(transcript.aa_change)
			output_file.write("\t")
			if gene.canonical_transcript_id==transcript.enstid: #Only print Impact prediction information for canonical transcript
				output_file.write(str(mutation.mutation_pred_rankscore_of_canonical))
			else:
				output_file.write("0")
			output_file.write("\t")
			if gene.canonical_transcript_id == transcript.enstid:#Only print Phastcons Conservation information for canonical transcript
				output_file.write(str(mutation.mutation_phastcons_of_canonical))
			else:
				output_file.write("0")
			output_file.write("\n")
		b=b+1

output_file.close()
