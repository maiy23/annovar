#!/usr/bin/perl
use strict;
use Getopt::Std;
use Data::Dumper;
use lib '/users/abcc/vuonghm/scripts.dir/perl_modules/JSON-2.90/lib/';
use JSON;
use FindBin;
use lib "$FindBin::Bin/";
use ABCC_utils;
use vars qw($opt_i $opt_n $opt_m $opt_o);
getopts("m:i:n:o:");
##for testing and laziness
my $dir="mydir";
$opt_m="$dir/mouse.coding.ProtPos.mapper";
$opt_i="$dir/20160624.novar";
$opt_n="failed.ProtPos";
my $newdir="newProvean";
##now it starts...
my $mapper_fn=$opt_m;#mapper file from original run 
my $new_ProtPos_fn=$opt_n;#new fn for prot pos for input back into provean_wrpr2.pl
$opt_o||="mm10";
my $taxon=getOrgType("$opt_o","taxonid");

open (MAPPER,"<$mapper_fn") or die "Cannot open $mapper_fn\n";
my %mapper;
while (<MAPPER>){
	my ($gene,$uniprot)=split("\t",$_);chomp $uniprot;
	$mapper{$uniprot}=$gene;
}
close MAPPER;
open (OUT,">$opt_n") or die "Cannot open $opt_n for writing\n";
open (FAILEDJOBS,"<$opt_i") or die "Cannot open $opt_i for reading\n";
while (<FAILEDJOBS>){
	my ($ldir,$uname)=split("/",$_);chomp $uname;$ldir=~s/ //;
	if (exists($mapper{$uname})){
		if (! -e "$dir/$ldir/$uname.var"){
			die "$dir/$ldir/$uname.var does not exist!\n";
		}
		my @res=split("\n",`cat $dir/$ldir/$uname.var`);
		foreach my $var (@res){
			print OUT $mapper{$uname} . ":$var\n";
		}
	}else{
		die "$uname does not exist!  Are you using the correct mapper file?...exiting\n";
	}
}
=head 
helper script for provean re-runs
After Provean_wrpr.pl fails for a large number of variants, run this script to map back the failed jobs list (file snippet below) back to 
[/bioinfoC/AVA/FDI/dbSNPtest/mouse]-> m mydir/20160624.novar 
 _A/A0A024QYR9
 _A/A0A087WSN1
 _A/A0A0A0MQK0
 _A/A0A0A6YWP8
 _A/A0A0G2JDF6
 ...
 Back to 
 [/bioinfoC/AVA/FDI/dbSNPtest/mouse]-> more mouse.coding.ProtPos                                                                                                         
Xkr4:p.D610D
Xkr4:p.P528P
Xkr4:p.G514G
Xkr4:p.F510F

To see which failed and rerun provean_wrpr2.pl with fasta aa pos check
=cut