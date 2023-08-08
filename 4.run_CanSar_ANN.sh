#!/bin/bash

# module load python/3.5.1

python CanSar_ANN_annotation_parser_v2.py \
--annotation_file big_example.vcf \
--maf_file TARGET.all.maf \
--id_mapping_file snpEff_genes.txt
