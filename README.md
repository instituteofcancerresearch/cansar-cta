# CanSAR - Cancer Target Association

## Summary
Code for calculating cancer target association scores, between a cancer type and target genes. Code has been re-written once and includes preprocessing steps to the VCF file. 

## Requirements

* SnpEff and SnpSift
* Python
* Perl (if preprocessing the vcf is needed)

## Processing stages (for reorganisation)
1. concatenate VCFs, write a 'maf' file to track patients to mutations and VAF, remove INFO field
2. run snpEff, sort and run snpSift
3. run CanSar_ANN_annotation_parser_v2.py
4. CanSAR_ANN_mutation_score_calculator_v2.py

## URL needed for last script 
http://www.ensembl.org/biomart/martview/1673cf0a9581742e160d9f9895153786?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id_version|hsapiens_gene_ensembl.default.feature_page.external_gene_name&FILTERS=hsapiens_gene_ensembl.default.filters.chromosome_name.%221,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y%22&VISIBLEPANEL=resultspanel

This URL should return mapping information from ENSG to ENST, a file needed for CanSAR_ANN_mutation_score_calculator_v2.py

## Bugra's original Readme

Input file with all mutations: mutations.vcf

### snpeff command
java -Xmx8g -jar snpEff.jar GRCh38.86 -v –canon mutations.vcf > mutations.annotated.onlycanonical.vcf

### sorting VCF lexicografically
cat mutations.annotated.onlycanonical.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "LC_ALL=C sort -k1,1 -k2,2n"}' > mutations.sorted.annotated.onlycanonical.vcf

### snpsift command
java -Xmx8g -jar SnpSift.jar dbnsfp -db ./snpEff_latest_core/snpEff/db/dbNSFP3.2a.txt.gz -f MetaSVM_score,MetaSVM_rankscore,phastCons100way_vertebrate -v mutations.sorted.annotated.onlycanonical.vcf > mutations.sorted.annotated.snpEffSnpSift.onlycanonical.vcf

### SCRIPT stage

### 1 - CanSar_ANN_annotation_parser.py

Input for the code: 
Id mapping file: snpeff_run.genes.txt (output of Snpeff run)
Mutation_file: mutations.sorted.annotated.snpEffSnpSift.onlycanonical.vcf
Maf_file: subset_of_maf_with_id_study_patient_genloc

Output:
annotation_summary_output_20191026.txt

### 2 - CanSAR_ANN_mutation_score_calculator.py
Input for the code:
Read_data: annotation_summary_output_20191026.txt
id_mapping_data: ENSG_to_ENST_to_GeneName.txt

Output of ranked list of scores:
Calculated_scores_for_CANSAR.txt

