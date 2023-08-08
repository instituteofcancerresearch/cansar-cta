#!/usr/bin/perl -w

use strict;

# Read a VCF file and split each line on tabs
# Set the array element [7] to '.' and write
# out the tab-joined fields to a new VCF

die "Error: First argument should be a path to a VCF file to strip the INFO string from\n" unless defined $ARGV[0];

my $in_vcf = $ARGV[0];
chomp $in_vcf;


my $out_vcf = $in_vcf;
if($out_vcf =~ /\.vcf$/){
	$out_vcf =~ s/vcf$/noInfo.vcf/;
}
else{
	$out_vcf .= '.noInfo.vcf';
}

open IN, "< $in_vcf" or die "Can't read $in_vcf: $!\n";
open OUT, "> $out_vcf" or die "Can't write to $out_vcf: $!\n";

while(<IN>){
	if(/^#/){
		print OUT;
		next;
	}
	my @fields = split(/\t/);
	$fields[7] = '.';
	print OUT join("\t", @fields);
}
close IN;
close OUT;


