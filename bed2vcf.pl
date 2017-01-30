#!/usr/bin/perl
use strict;
my $input=$ARGV[0];
my $output=$ARGV[1];
$output||="$input.vcf";
my $headers;
getHeaders();

open (INPUT,"<$input") or die "Cannot open $input\n";
open (OUT,">$output") or die "Cannot write to $output\n";
open (TMP,">tmp") or die "Canno open tmp\n";
print OUT $headers;
close OUT;
my $line=1;
while (<INPUT>){
	#chr,start,stop,$ref,$var...
	my ($chr,$start,$stop,$ref,$var,@rest)=split("\t",$_);
	chomp $var;
	if ($#rest>-1){
		chomp $rest[$#rest] ;
	}else{
		$rest[0]=".";
	}
	$chr=~s/chr//;
	if ($start==$stop && $ref=~/^[ATCG\-]$/ && $var=~/^[ATCG]$/){
		print TMP join ("\t",$chr,$start,".",$ref,$var,".",".",".",join("|",@rest,"ABCC_LN=$line"))."\n";
	}else{
		print TMP join ("\t",$chr,$start,".",$ref,$var,".",".",".",join("|",@rest,"ABCC_LN=$line"))."\n";
	}
	$line++;

}
close TMP;
system ("sort -k1,1n -k2,2n tmp |grep -ve '^#' >>$output\n");

sub getHeaders(){
	$headers="##fileformat=VCFv4.0\n".
	"##AVIAv2.0\n".
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}