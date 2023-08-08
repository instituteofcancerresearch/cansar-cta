

module load java/sun8/1.8.0u66
module load snpeff/4.3t

DIR='./'
VCF_PRFX='TARGET.all'

perl strip_info_from_vcf.pl ${DIR}${VCF_PRFX}.vcf 

java -Xmx8g \
-jar /apps/snpeff/4.3t/snpEff/snpEff.jar \
-v -canon GRCh38.86 \
${DIR}${VCF_PRFX}.noInfo.vcf  > ${DIR}${VCF_PRFX}.noInfo.annotated.onlycanonical.vcf

cat ${DIR}${VCF_PRFX}.noInfo.annotated.onlycanonical.vcf \
| awk '$1 ~ /^#/ {print $0;next} {print $0 | "LC_ALL=C sort -k1,1 -k2,2n"}' \
> ${DIR}${VCF_PRFX}.noInfo.annotated.onlycanonical.sorted.vcf

java -Xmx8g -jar /apps/snpeff/4.3t/snpEff/SnpSift.jar dbnsfp \
-db ${DIR}dbNSFP4.0a.txt.gz \
-f MetaSVM_score,MetaSVM_rankscore,phastCons100way_vertebrate \
-v ${DIR}${VCF_PRFX}.noInfo.annotated.onlycanonical.sorted.vcf \
> ${DIR}${VCF_PRFX}.noInfo.annotated.onlycanonical.sorted.snpSift.vcf

