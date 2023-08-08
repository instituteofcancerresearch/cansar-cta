#!/usr/bin/perl -w

use strict;

# script to retrofit the VCF file for Bugra's Mutation Analysis
# script with the variant allele fraction (alt reads / total depth)
# in column 3 (ID). The value between 0 and 1 should be a suffix to
# the mutaton index value, separated by an underscore
# e.g. 1_0.2

# VCFs have the following columns - need to get format from FORMAT [8] and AD from TUMOUR [10]
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR

my $vcf_in = 'TARGET.all.noInfo.annotated.onlycanonical.sorted.snpSift.vcf';

my $tum_samp_col_idx = 10;

my $vcf_out = $vcf_in;
if($vcf_out =~ /\.vcf$/){
	$vcf_out =~ s/vcf/addedVaf.vcf/;
}
else{
	$vcf_out .= '.addedVaf.vcf';
}

open IN, "< $vcf_in" or die "Can't read $vcf_in: $!\n";
open OUT, "> $vcf_out" or die "Can't write to $vcf_out: $!\n";

while(<IN>){
	if(/^#/){print OUT; next;}
	# Assuming here that col index 10 is the tumour
	# and the second :-sep field is AD=ref,alt
	my @fields = split(/\t/);
	my @samp = split(/:/, $fields[$tum_samp_col_idx]); # col idx 10 == col 11
	unless($samp[1] =~ /\d+,\d+/){
		die "unable to understand the AD field $samp[1] - expected two ints separated by a comma\n";
	}
	my ($ref, $alt) = split(/,/, $samp[1]);
	my $vaf = $alt / ($alt + $ref);
	$fields[2] = $fields[2] . "_" . $vaf;
	print OUT join("\t", @fields);
}

close IN;
close OUT;




