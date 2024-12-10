#!/usr/bin/perl

# # # # # #
# filterVCFonQuality.pl
# script obtained from Krysia Nadachowska-Brzyska in 2015 

# ====================================================
# Takes an VCF file and masks bases with N if they
# don't fulfill certain quality and depth criterias. 
# ====================================================


# Example 
# samtools mpileup -uf file.fasta sample.bam |bcftools view -cg - | perl filterVCFonQuality.pl /dev/stdin 25 20 19.8 >new.vcf

use strict;
use warnings;

# Input parameters
my $VCF = $ARGV[0];		# Could be a real vcf file, or piped output from samtools/bcftools
my $MQ = $ARGV[1];		# mapping quality threshold	#Normally set to 25
my $FQ = $ARGV[2];		# consensus quality threshold	# Normally set to 20
my $AVGDEPTH = $ARGV[3];	# average depth for the sample (thresholds used: avgd/3< && <avgd*2)


# Go through the vcf file
open(FILE, $VCF);
my ($cnt, $lm, $lb, $ld, $hd)  = (0,0,0,0,0);
while(<FILE>) {

	if(/^#/) {	# Comment lines are printed
		print $_; 
	}
	else {
		# Create an array with all columns
		my @arr = split(/\s+/, $_);

		# Set initial quality to zero (updated below if they exist in the file)
		my ($dp, $mq, $fq) = (0, 0, 0);

		if ($arr[7] =~ m/DP\=(\d+)/){
			$dp=$1;
		}
		if ($arr[7] =~ m/MQ\=(\d+)/){
			$mq=$1;
		}
		if ($arr[7] =~ m/FQ\=(\d+)/){
			$fq=$1;
		}
		if ($arr[7] =~ m/FQ\=(-\d+)/){
			$fq=$1;
		}
		# Check if thresholds are fulfilled (and save stats)
		if ($mq < $MQ){
			$lm++;
			$arr[3]="N";
			$arr[4]=".";
		}
		#if($dp < $AVGDEPTH/3) {
		if($dp < 10) {
			$ld++;
			$arr[3]="N";
			$arr[4]=".";
		}
		if($dp>$AVGDEPTH*2) {
			$hd++;
			$arr[3]="N";
			$arr[4]=".";
		}
		if (abs($fq) < $FQ){
			$lb++;
			$arr[3]="N";
			$arr[4]=".";
		} 

		print join("\t", @arr), "\n";
		$cnt++;
	}
}
close(FILE);

# Print stats
print STDERR "Out of $cnt VCF rows, $lm had too low mapping quality, $lb had too low consensus quality, $ld had too low and $hd had too high depth respectively.\n";


