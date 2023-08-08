#!/usr/bin/perl -w

use strict;

# Concatenate VCFs for Bugra's scripts...
# Read file location and IDs from manifest
# Open each VCF file and write to a new file
# skip lines begining with #
# skip lines without PASS in filter column
# increment a counter snf store it as an index in col3 of catted VCF
# add index, study, patient and gcoords to 'maf' file

my $sample_sheet_file = 'gdc_sample_sheet.2020-03-16.tsv';
my $gdc_dir = 'gdc_download_20200316_093641.459183/';
my $output_vcf_file = 'TARGET-all.vcf';
my $output_maf_file = 'TARGET-all.maf';


open VCFOUT, "> $output_vcf_file" or die "Can't write to output VCF $output_vcf_file: $!\n";
open MAFOUT, "> $output_maf_file" or die "Can't write to output MAF $output_maf_file: $!\n";

open SS, "< $sample_sheet_file" or die "Can't read sample information from $sample_sheet_file: $!\n";
my $header= <SS>;

my $vcf_row_index = 0;
while(<SS>){
	my @fields = split(/\t/);
	my $dir = $fields[0];
	my $vcffile = $fields[1];
	$vcffile =~ s/\.gz$//;
	my $project = $fields[4];
	my $patient = $fields[5];
	$patient =~ s/,.*//;
	
	open IN, "< ${gdc_dir}/${dir}/$vcffile" or die "Can't read ${dir}/$vcffile: $!\n";
	while(<IN>){
		next if /^#/;
		my @fields = split(/\t/);
#		print "Looking at $fields[6]\n";
		next unless $fields[6] eq "PASS";
		$vcf_row_index ++;
		print "$vcf_row_index\n";
		$fields[2] = $vcf_row_index;
		print VCFOUT join("\t", @fields);
		my $position = "$fields[0]" . ":" . "$fields[1]";
		print MAFOUT "$vcf_row_index $project $patient $position\n";
	}
	
}

close SS;
close VCFOUT;
close MAFOUT;



