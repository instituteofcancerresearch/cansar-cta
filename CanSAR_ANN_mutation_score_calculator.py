import os
import statistics
import pandas as pd
from collections import Counter

class Gene:
    def __init__(self):
        self.genename=""
        self.ensemblg_id=""
        self.ensemblt_id=""
        self.mutation_types=[]
        self.mutation_pred_rankscore=[]
        self.mutation_phastcons=[]
        self.mutation_impact=[]
        self.mutation_depth=[]
        self.mutation_chr_pos = []
        self.mutation_protein_pos=[]
        self.mutation_studies=[]
        self.mutation_patients=[]
        self.most_freq_mut_percentage_in_cohort=[]
        self.percentage_of_cohort_with_mut_inthisgene=[]
        self.mutation_impact_list=[]
        self.median_impact_score=0
        self.mutation_overall_score=[]


#READ THE METADATA produced using CanSar_ANN_annotation_parser.py
file='annotation_summary_output_20191026.txt'

#Use id mapping information from ENSG TO ENST and create dictionary
id_mapping_data='ENSG_to_ENST_to_GeneName.txt'
idmapping_read=pd.read_csv(id_mapping_data,sep="\t")

idmapping_read=idmapping_read.loc[:,['Gene stable ID','Transcript stable ID','Gene name']]
enst_to_ensg_dict={}
genename_to_ensg_dict={}
ensg_to_genename_dict={}
for i in idmapping_read.index:
    enst_to_ensg_dict[idmapping_read.at[i,'Transcript stable ID']]=idmapping_read.at[i,'Gene stable ID']
    genename_to_ensg_dict[idmapping_read.at[i,'Gene name']]=idmapping_read.at[i,'Gene stable ID']
    ensg_to_genename_dict[idmapping_read.at[i,'Gene stable ID']]=idmapping_read.at[i,'Gene name']
#END of dictionary creation for ID mapping

file_read=open(file,"r")

lines=file_read.readlines()

annotated_genes=[]
linecount=0
studylist=[]

#START PARSING THE METAFILE
for line in lines:
    if line.startswith("EnsembGenomeID")==False: #ignore the header
    ###if line.startswith("Gene stable ID")==False: #ignore the header
        columns=line.split("\t")
        newgene = Gene()
        newgene.ensemblg_id=columns[0]
        newgene.genename = columns[1]
        newgene.ensemblt_id = columns[2]
        # ASSIGN ENST, ENSG, Genename info properly
        if newgene.ensemblt_id in enst_to_ensg_dict.keys():
            newgene.ensemblg_id = enst_to_ensg_dict[newgene.ensemblt_id]
            newgene.genename=ensg_to_genename_dict[newgene.ensemblg_id]
        elif newgene.ensemblg_id in genename_to_ensg_dict.keys():
            newgene.ensemblg_id = genename_to_ensg_dict[newgene.genename]

        numbernottoadd=[]
        index_2=0
        for string in columns[8].split(","):
            if string=="synonymous_variant" or string=="TF_binding_site_variant" or string=="protein_protein_contact" or string=="TFBS_ablation" or string== "intron_variant": # INTRONIC or unreliable MUTATION TYPES
                numbernottoadd.append(index_2) #if mutation type is one these, then add the index to DONOTADD list.

            index_2=index_2+1

        
        notadd = 0
        for string in columns[9].split(","):#MUTATION PROTEIN POSITIONS
            if string!="" and notadd not in numbernottoadd :
                newgene.mutation_protein_pos.append(int(string))
            elif notadd not in numbernottoadd:
                numbernottoadd.append(notadd)
            notadd=notadd+1

        notadd = 0
        for string in columns[3].split(","):#PATIENTS_W_MUTATIONS
            if string!="" and notadd not in numbernottoadd :
                newgene.mutation_patients.append(string)
            elif notadd not in numbernottoadd:
                numbernottoadd.append(notadd)
            notadd=notadd+1
        notadd = 0
        for string in columns[4].split(","):#STUDIES_W_MUTATIONS
            if string!="" and notadd not in numbernottoadd :
                newgene.mutation_studies.append(string)
            elif notadd not in numbernottoadd:
                numbernottoadd.append(notadd)
            notadd=notadd+1

            if string not in studylist and string != "":
                studylist.append(string)
        notadd = 0
        for string in columns[5].split(","):#MUTATION IMPACT PREDICTION SCORE
            if string!="" and notadd not in numbernottoadd :
                if string==".":
                    newgene.mutation_pred_rankscore.append(0)
                else:
                    newgene.mutation_pred_rankscore.append(float(string))
            elif notadd not in numbernottoadd:
                numbernottoadd.append(notadd)
            notadd=notadd+1

        notadd = 0
        for string in columns[6].split(","):#PHASTCONS CONSERVATION SCORES
            if string!="" and notadd not in numbernottoadd :
                newgene.mutation_phastcons.append(float(string))
            elif notadd not in numbernottoadd:
                numbernottoadd.append(notadd)
            notadd=notadd+1

        notadd = 0
        for string in columns[7].split(","):#MUTATION CHROMOSOME LOCATIONS
            if string!="" and notadd not in numbernottoadd :
                newgene.mutation_chr_pos.append(string)
            elif notadd not in numbernottoadd:
                numbernottoadd.append(notadd)
            notadd=notadd+1

        notadd = 0
        for string in columns[8].split(","):#MUTATION TYPES
            if string!="" and notadd not in numbernottoadd :
                newgene.mutation_types.append(string)
            elif notadd not in numbernottoadd:
                numbernottoadd.append(notadd)
            notadd=notadd+1

        notadd = 0
        for string in columns[10].rstrip().split(","):#MUTATION IMPACTS (LOW,MODERATE,HIGH)
            if string!="" and notadd not in numbernottoadd :
                newgene.mutation_impact.append(string)
            elif notadd not in numbernottoadd:
                numbernottoadd.append(notadd)
            notadd=notadd+1


        notadd = 0
        for string in columns[11].split(","): #MUTATION DEPTH
            if string!="" and notadd not in numbernottoadd :
                newgene.mutation_depth.append(float(string))
            elif notadd not in numbernottoadd:
                numbernottoadd.append(notadd)
            notadd=notadd+1


        if len(newgene.mutation_depth)>0: # CALCULATE THE 4th COMPONENT OF THE SCORE, integrating CONSERVATION SCORE,IMPACT PREDICTION and MUTATION DEPTH
            for a in range(0,len(newgene.mutation_studies)):
                if newgene.mutation_impact[a]=="HIGH":
                    newgene.mutation_impact_list.append((newgene.mutation_phastcons[a]/5)+(newgene.mutation_pred_rankscore[a]/5)+(newgene.mutation_depth[a]/5)+0.4)
                elif newgene.mutation_impact[a]=="MODERATE":
                    newgene.mutation_impact_list.append((newgene.mutation_phastcons[a] / 5) + (newgene.mutation_pred_rankscore[a] / 5) + (newgene.mutation_depth[a] / 5) + 0.3)
                elif newgene.mutation_impact[a]=="LOW":
                    newgene.mutation_impact_list.append((newgene.mutation_phastcons[a] / 5) + (newgene.mutation_pred_rankscore[a] / 5) + (newgene.mutation_depth[a] / 5) + 0.15)
                else:
                    newgene.mutation_impact_list.append((newgene.mutation_phastcons[a] / 5) + (newgene.mutation_pred_rankscore[a] / 5) + (newgene.mutation_depth[a] / 5) + 0.10)

            newgene.median_impact_score=statistics.median(newgene.mutation_impact_list) # Get the Median score for this GENE's impact scores.

            annotated_genes.append(newgene)

    print(linecount)
    linecount=linecount+1

#ALL MUTATION INFORMATION FROM METADATA ARE STORED IN ANNOTATED_GENES LIST

study_burden=[]
study_percentage_of_recurrent_mutations_list=[]
#START OF 2nd COMPONENT ( RECURRENCY) and 3rd COMPONENT CALCULATION (MUTATION BURDEN)
for study in studylist: # for the studies in the vcf file
    print(study)
    study_mutationlist=0
    study_patientlist=[] # initialize patients in this study
    no_of_mutated_genes_in_that_study=0 # initialize no of mutated genes in this study
    no_of_recurrent_mutations_in_that_study=0 # initialize number of recurrent mutations in this study
    for gene in annotated_genes:
        for a in range(0,len(gene.mutation_studies)):
            if gene.mutation_studies[a]==study:
                study_mutationlist=study_mutationlist+1 # no.of mutations in this study
            if gene.mutation_studies[a]==study and gene.mutation_patients[a] not in study_patientlist:
                study_patientlist.append(gene.mutation_patients[a]) # Add patients to the list if they haven't added already. This should ultimately give you the list of unique patients which have mutations in this study type

    study_burden.append(study_mutationlist/len(study_patientlist)) # study burden = no.of mutations in that study / no. of unique patients with mutation in that study
    for gene in annotated_genes: #for all genes
        patients_coming_from_thisgene = [] #initialize patients with mutation in this gene
        mutation_positions=[] #initialize mutation positions for recurrency calculations
        check=False
        for a in range(0,len(gene.mutation_studies)): #for all studies in TCGA
            if gene.mutation_studies[a]==study and gene.mutation_patients[a] not in patients_coming_from_thisgene: #Count unique number of patients in this study that has mutation in this specific gene
                patients_coming_from_thisgene.append(gene.mutation_patients[a])
            if gene.mutation_studies[a]==study:
                if gene.mutation_impact[a]=="HIGH" and max(gene.mutation_protein_pos)<7500: #IF the mutation is Frameshift,Stoploss,Stopgain and protein is not longer than 7500, then add the position of mutation as "position 0". This is to avoid different frameshift mutation not being counted as recurrent. 7500 limit is to remove the artefacts like TTN.
                    mutation_positions.append(0)
                else:
                    mutation_positions.append(gene.mutation_protein_pos[a])
                check=True

        if check==True: #For recurrency calculation
            result = most_common, num_most_common = Counter(mutation_positions).most_common(1)[0] #calculate most recurrent position

            if result[1]>=1: # if there is at least 1 position with 2 mutations
                no_of_recurrent_mutations_in_that_study=no_of_recurrent_mutations_in_that_study+1
                no_of_mutated_genes_in_that_study = no_of_mutated_genes_in_that_study + 1
                if result[1]>=40:#Give higher weight to recurrence >40
                    gene.most_freq_mut_percentage_in_cohort.append(2)
                else: #Give importance in relation with 20 mutations
                    gene.most_freq_mut_percentage_in_cohort.append(result[1]/20)# / len(mutation_positions))
            else: #If there isn't any recurrent mutation, give less importance
                no_of_mutated_genes_in_that_study = no_of_mutated_genes_in_that_study+1
                gene.most_freq_mut_percentage_in_cohort.append(0.05)
        else:
            gene.most_freq_mut_percentage_in_cohort.append(0) # if there isn't any mutation then assign "0" score
        gene.percentage_of_cohort_with_mut_inthisgene.append(len(patients_coming_from_thisgene) / len(study_patientlist)) #Patients with mutation in this gene divided by Number of Unique Patients in that study

    study_percentage_of_recurrent_mutations_list.append(no_of_recurrent_mutations_in_that_study / no_of_mutated_genes_in_that_study) # Percentage of recurrent mutation positions to total number of mutation positions

study_burden_percentages=[]
study_percentage_of_recurrent_mutations_percentages=[]
for value in study_percentage_of_recurrent_mutations_list: #Convert per study recurrent mutation percentages to relative amount
    study_percentage_of_recurrent_mutations_percentages.append(value/sum(study_percentage_of_recurrent_mutations_list))

for value in study_burden:
    study_burden_percentages.append(value/sum(study_burden)) #Convert per study recurrent mutation percentages to relative amount

normalised_study_burden_percentages=[]
normalised_study_percentage_of_recurrent_mutations=[]

for value in study_percentage_of_recurrent_mutations_percentages:#Normalize recurrent mutation percentages - numbers are in accordance with the final computation score
    if value>statistics.mean(study_percentage_of_recurrent_mutations_percentages)+statistics.stdev(study_percentage_of_recurrent_mutations_percentages): # mean + 1 standard deviation
        normalised_study_percentage_of_recurrent_mutations.append(0.2125)
    elif value>statistics.mean(study_percentage_of_recurrent_mutations_percentages)+(2*statistics.stdev(study_percentage_of_recurrent_mutations_percentages)):# mean + 2 standard deviation
        normalised_study_percentage_of_recurrent_mutations.append(0.425)
    else:
        normalised_study_percentage_of_recurrent_mutations.append(0.05) #else assign low score

for value in study_burden_percentages:#Normalize the study burden percentages - numbers are in accordance with the final computation score
    if value > statistics.mean(study_burden_percentages) + statistics.stdev(study_burden_percentages): # mean + 1 standard deviation
        normalised_study_burden_percentages.append(0.2125)
    elif value > statistics.mean(study_burden_percentages) + (2 * statistics.stdev(study_burden_percentages)): # mean + 2 standard deviation
        normalised_study_burden_percentages.append(0.425)
    else:
        normalised_study_burden_percentages.append(0.05)  #else assign low score


#START OF OUTPUTTING THE FILE

output_file=open("Calculated_scores_for_CANSAR.txt","w")

#Header Line
output_file.write("EnsemblID"+"\t"+"GeneName"+"\t"+"StudyName"+"\t"+"Gene Mutation Score"+"\t"+"% of cohort with mutation in this gene"+"\t"+"Most Frequently mutated position vs Total number of mutations"+"\t"+"Cohort Mutational Burden"+"\t"+"Average number of mutation per patient in that cohort"+"\t"+"% of mutated genes with recurrent mutations"+"\t"+"Gene Median Impact Score")
output_file.write("\n")

counter=0
#Print the association score information
for gene in annotated_genes:
    #if gene.ensembl_id=="ENST00000254227":
        #letmeknow=True
    print(counter)
    for ind in range(0, len(studylist)):
        gene.mutation_overall_score.append(
            (gene.percentage_of_cohort_with_mut_inthisgene[ind]+1) * (gene.most_freq_mut_percentage_in_cohort[ind]) * (1.85 - (
                        normalised_study_burden_percentages[ind] + normalised_study_percentage_of_recurrent_mutations[
                    ind])) * (gene.median_impact_score+1))
        #gene.mutation_overall_score.append(gene.percentage_of_cohort_with_mut_inthisgene[ind]*gene.most_freq_mut_percentage_in_cohort[ind]*(1-(normalised_study_burden_percentages[ind]+normalised_study_percentage_of_recurrent_mutations[ind]))*gene.median_impact_score)
        output_file.write(gene.ensemblt_id+"\t"+gene.genename+"\t"+studylist[ind]+"\t"+str(gene.mutation_overall_score[-1])+"\t"+str(1+gene.percentage_of_cohort_with_mut_inthisgene[ind])+"\t"+str(gene.most_freq_mut_percentage_in_cohort[ind])+"\t"+str(1.85-(normalised_study_burden_percentages[ind]+normalised_study_percentage_of_recurrent_mutations[ind]))+"\t"+str(normalised_study_burden_percentages[ind])+"\t"+str(normalised_study_percentage_of_recurrent_mutations[ind])+"\t"+str(1+gene.median_impact_score))
        output_file.write("\n")

    counter=counter+1
output_file.close()

