#!/usr/bin/perl
use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use ABCC_utils;
use vars qw($opt_i $opt_f $opt_a $opt_d $opt_D);
getopts("f:i:d:D:a:");
my $usage=qq(
	$0 
	[REQUIRED]
		-i <bed input file name>
		-a <folder Name> 
			Outputs results into this directory; makes directory if it does not exist
	[OPTIONAL]
		-d <database folder> if different than mouse
		-D <UCSC build> default mm10;
);
my $db_folder="/SeqIdx/annovardb/SIFT4G-Annotator/Database/GRCm38.74/";
# my $sift_output_file=
if (defined $opt_d){
	$db_folder=$opt_d;
}
if (!defined $opt_D){
	$opt_D="mm10";
}
if (!defined $opt_a || !defined $opt_i){
	print "You must supply in input file and directory\n\n$usage\n";exit;
}

my $bin=`dirname $0`;chomp $bin;
my $bed2vcf_exe="$bin/bed2vcf.pl";
my $perl=locate("perl");
my $sift4g_annotator_exe="java -Xmx512m -jar /SeqIdx/annovardb/SIFT4G-Annotator/SIFT4G_Annotator_v2.4.jar -c -i $opt_i.vcf -d  $db_folder -t -r $opt_a";

#first take input and convert back to VCF file
#SIFT4G requires a VCF file
###REINSTATE
my $success=SAFE_PROCESS("$perl $bed2vcf_exe $opt_i $opt_i.vcf\n");
if (!$success){
	die "Could not execute $bed2vcf_exe at ". __LINE__." in $0\n";
}
print "about to run $sift4g_annotator_exe at line " .__LINE__."\n";
$success=SAFE_PROCESS("$sift4g_annotator_exe");
if (!$success){
	die "Could not run $sift4g_annotator_exe at line ". __LINE__. " in $0\n";
}

my $grep=locate("grep");
my $cut=locate("cut");
rearrange();
unlink ("$opt_a/$opt_i\_SIFTpredictions.vcf", "$opt_a/$opt_i\_SIFTannotations.xls","tmp");#remove intermediate files

sub rearrange{
	open (OUT,">$opt_i.$opt_D\_sift4G\n") or die "Cannot open $opt_i.$opt_D\_sift4G\n";
	print OUT "#sift4G\n";
	print "In rearrange\n";
	open (IN,"<$opt_a/$opt_i" . "_SIFTannotations.xls") or die "Cannot read $opt_a/$opt_i" . "_SIFTannotations.xls\n";
	my $sorted_line_num=0;my %mapper;
	while (my $line=<IN>){
		$line=~s/[\n\r]$//;chomp $line;
		next if ($line!~/^[\dXYM]/ || $line=~/^\s*$/);
		my @lines=split("\t",$line);
		my $loc=join(",",@lines[0..3]);
		$mapper{$loc}=$lines[16].";".$lines[12];
	}
	close IN;
	open (ORIG,"<$opt_i") or die "Cannot open $opt_i\n";
	while (my $line=<ORIG>){
		next if ($line=~/^#/ || !$line);
		my @lines=split("\t",$line);
		my $loc=join(",",@lines[0,1,3,4]);
		if (exists ($mapper{$loc})){
			print OUT "$mapper{$loc}\n";
		}else{
			print OUT "-\n";
		}
	}
	close ORIG;close OUT;
}
