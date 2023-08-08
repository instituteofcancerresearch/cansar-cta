import os
import bisect
import argparse

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
#     the mutation index number, seq-project, patient-id and mutation location.


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





#Location of file to map ENS ids with each other
#id_mapping_file="/Users/bozer/Documents/snpEff_latest_core/snpEff/20191018_mac_home_run.genes.txt" #"/Users/bozer/Documents/SEQUENCING/snpEff_latest_core_WORKINGLATESTONE/snpEff/TCGA_snpeff_canonical_onlyprotein.genes.txt"#"/Users/bozer/Documents/snpEff_latest_core/snpEff/2019101_mac_home_run.genes.txt"#

file_id_mapping=open(id_mapping_file,"r")
id_dict={}
#create dictionary of IDs, from ENST to ENSG
for idline in file_id_mapping:
    if idline.startswith("#")==False:
        columns=idline.split("\t")
        id_dict[columns[1]]=columns[2]

#Location of dbNSFP-snpeffed annotated mutation file with all TCGA mutations
#annotation_file="/Volumes/Samsung_T5/CANSAR/TCGA_mutations_SnpEffed_20191014/TCGA_mutations_withIDmapping_onlyprotein_sorted_dbNSFPadded_20191014.vcf"
#Location of MAF file which has been used as identifier
#maf_file="/Volumes/Samsung_T5/CANSAR/subset_of_maf_with_id_study_patient_genloc.txt"

file_maf=open(maf_file,"r")
file_annotation=open(annotation_file,"r")

maf_lines=file_maf.readlines()

#annotated_genes=[]
ensemblgnamelist=[]
#true_genenamelist=[]
genelist=[]

linenumber=0
#Do not consider these mutations as these are either non-coding or not well defined
unwanted_mutation_types=["TF_binding_site_variant","protein_protein_contact","TFBS_ablation","intron_variant","intragenic_variant","structural_interaction_variant","non_coding_transcript_exon_variant"]


#READING THE MUTATION VCF FILE and LINKING TO MAF FILE INFORMATION
with open(annotation_file,"r") as f:
    for line in f:
        if line.startswith("#") == False: #read each line which doesn't start with # (header starts with #)
            line=line.rstrip("\n")
            if line.split("\t")[1]!="0": #if chromosome number is not "0", do the the below
                print(linenumber)
                #initialize some of the variants that will change with every line
                mutation_rankscore=0
                mutation_phastcons=0
                line_parser=line.split(";")
                aachangesingleletter = ""

                newmutation=Mutation()

                for string in line.split(";"): # Parse annotation field of the line
                    if string.startswith("dbNSFP_MetaSVM_rankscore"): #Rankscore=Mutation Impact Prediction
                        temp_score=string.split("=")[1].split(",")
                        for score in temp_score:
                            if score!=".":
                                mutation_rankscore=score#string.split("=")[1]
                                break
                            else:
                                letmeknow=True
                    elif string.startswith("dbNSFP_phastCons100way_vertebrate"): # Get PhastCons Conservation Score
                        temp_phastcons=string.split("=")[1].split(",")
                        for phastcons in temp_phastcons:
                            if phastcons!=".":
                                mutation_phastcons=phastcons#string.split("=")[1]
                                break
                            else:
                                letmeknow=True

                    elif string.startswith("dbNSFP_aaref"): #Reference Aminoacid
                        aachangesingleletter=aachangesingleletter+string.split("=")[1]
                    elif string.startswith("dbNSFP_aapos"): # Aminoacid Position of Mutation
                        aachangesingleletter=aachangesingleletter+string.split("=")[1]
                    elif string.startswith("dbNSFP_aaalt"): #Aminoacid result of Mutation
                        aachangesingleletter=aachangesingleletter+string.split("=")[1]

                newmutation.aa_change_singleletter=aachangesingleletter
                newmutation.mutation_pred_rankscore_of_canonical = mutation_rankscore
                newmutation.mutation_phastcons_of_canonical = mutation_phastcons

                #Parse Different Transcript information for same Mutation
                transcripts=line.split(";")[0].split(",")
                
                print("length of transcripts: %d" % (len(transcripts)))
                
                # ? Don't get what this shoul do.
                # transcripts contains 7 elements.
                if len(line.split(";"))>1:
                    if len(transcripts)>1:
                        letmeknow=True
                if len(transcripts)>0:
                    for transcript in transcripts:

                        columns = transcript.split("\t")
                        if len(columns)>1: # if there is annotation information for that transcript
                            info_field = columns[7].split(";") #parse annotation field
                            mutation_id = columns[2].split("_")[0] #mutation ID that links to MAF file
                            mutation_depth = float(columns[2].split("_")[1]) #Retrieve Mutation Depth
                            
                            newmutation.mutation_depth=mutation_depth
                            mutation_index = int(mutation_id)
                            if mutation_index != int(maf_lines[mutation_index].split(" ")[0]):
                                letmeknow = True

                            study_id = maf_lines[mutation_index].split(" ")[1] #Get study information from Maf file
                            patient_id = maf_lines[mutation_index].split(" ")[2].rstrip() #TCGA patient ID
                            mutation_chr_pos = maf_lines[mutation_index].split(" ")[3].rstrip() #Chromosome position of mutation
                        else: # if there isn't detailed annotation information just use the first column information
                            info_field=columns[0].split(";")

                        if info_field[0].split("|")[1] not in unwanted_mutation_types: # if mutation type is among the ones we are interested in
                            newtranscript = Transcript()

                            mutation_protein_pos = ""
                        #for i in range(1, len(info_field[0].split("|")), 15):
                        #if info_field[0].split("|")[1] not in unwanted_mutation_types: #Among the annotations, choses the one that has biggest impact and not in unwanted mutation type list
                            
                            
                            #parse the annotation field coming from SNPEFF
                            impact = info_field[0].split("|")[2]
                            ensembl_id = info_field[0].split("|")[4].split(".")[0] #ENSG id
                            enst_id=info_field[0].split("|")[6].split(".")[0] #ENST id
                            mutation_type = info_field[0].split("|")[1]
                            genename = info_field[0].split("|")[3]

                            #DECIDE ON POSITION

                            pos_to_add = info_field[0].split("|")[12].split("/")[0]

                            mutation_protein_pos = pos_to_add

                            if impact!="MODIFIER": #Create new transcript instance if impact is LOW,HIGH or MODERATE as MODIFIER ones are in intronic-intergenic regions
                                newtranscript.enstid = enst_id if enst_id.startswith("ENST") else ""
                                newtranscript.mutation_type = mutation_type
                                newtranscript.mutation_impact = impact
                                newtranscript.mutation_protein_pos = mutation_protein_pos
                                newtranscript.aa_change = info_field[0].split("|")[10]

                                #add mutation information of this transcript to the mutation class
                                newmutation.transcripts.append(newtranscript)
                                newmutation.mutation_chr_pos = mutation_chr_pos

                                if ensembl_id not in ensemblgnamelist:#if mutation in this Gene is not added before create new Gene instance
                                    newgene = Gene()
                                    newgene.genename = genename
                                    newgene.ensg_id = ensembl_id
                                    newgene.mutations.append(newmutation)
                                    newgene.mutation_patient.append(patient_id)
                                    newgene.mutation_study.append(study_id)
                                    newgene.mutposlist.append(mutation_chr_pos)

                                    if ensembl_id in id_dict.keys():
                                        newgene.canonical_transcript_id=id_dict[ensembl_id] #find canonical transcipt through the dictionary file that was defined in the beginning
                                    else:
                                        print(ensembl_id) #if ensembl ID is not found then print those IDs.

                                    genelist.append(newgene) #add the newgene to List of Genes
                                    ensemblgnamelist.append(ensembl_id) # add Gene name to List of Genenames for indexing purposes
                                else:#if there is already mutation in this gene

                                    genename_index = ensemblgnamelist.index(ensembl_id) #find index of the gene so that new mutation can be added to Gene
                                    check=False

                                    if mutation_chr_pos not in genelist[genename_index].mutposlist: # if mutation position is already added before then Check=True
                                        check=True
                                    else: # if the new mutation location is not already added before
                                        #posindex=genelist[genename_index].mutposlist.index(mutation_chr_pos)
                                        posindex = len(genelist[genename_index].mutposlist) - genelist[genename_index].mutposlist[::-1].index(mutation_chr_pos) - 1
                                        if patient_id!=genelist[genename_index].mutation_patient[posindex]: #if mutation position is already there but if mutation belongs to another patient then make Check=True
                                            check=True
                                        else:
                                            check=False

                                    if check==True: #if this mutation is not added before per checks in previous stage
                                        genelist[genename_index].mutations.append(newmutation)
                                        genelist[genename_index].mutation_patient.append(patient_id)
                                        genelist[genename_index].mutation_study.append(study_id)
                                        genelist[genename_index].mutposlist.append(mutation_chr_pos)

            else:
                print("Faulty Line at:"+str(linenumber)) #Print the line number if line is not as expected

        linenumber = linenumber + 1


#PARSING FINISHED. BELOW IS FOR CREATING OUTPUT FILES

output_file=open("annotation_summary_output_20191026.txt","w") #FILE NAME FOR OUTPUT

#PRINT HEADER FIRST
output_file.write("EnsembGenomeID"+"\t"+"Gene Name"+"\t"+"EnsemblTranscriptID"+"\t"+"Mutations"+"\t"+"Patients"+"\t"+"Studies"+"\t"+"Pred_Rankscore"+"\t"+"PhastCons"+"\t"+"Chr Pos"+"\t"+"Mutation Type"+"\t"+"Protein Pos"+"\t"+"Mutation Impact"+"\t"+"Mutation Depth")
output_file.write("\n")

#PRINT THE METAFILE FOR EACH GENE
for gene in genelist:

    mutnottoadd=[]
    mutindex = 0
    for mut in gene.mutations:
        addtooutput = False
        for enst in mut.transcripts:
            if enst.enstid == gene.canonical_transcript_id:
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

    b = 0
    for patient in gene.mutation_patient:
        if b not in mutnottoadd:
            output_file.write(patient)
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
