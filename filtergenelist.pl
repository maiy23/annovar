#!/usr/bin/perl
use strict;
use FindBin;
use Getopt::Std;
use Data::Dumper;
use FindBin;
use Tie::File;
use lib "$FindBin::Bin/../perl_modules/Config-IniFiles-2.38";#config
use lib "$FindBin::Bin/";
use Config::IniFiles;
use lib "/bioinfoC/hue/annovar/annovar_dev/Parallel-ForkManager-0.7.9/lib";
use lib "/bioinfoC/AVA/prod_scripts/util/perl_modules";#parse and sendmail
use parse;
use vars qw($opt_c $opt_f $opt_i);
 getopts("f:i:c:");
 umask(0022);
 die qq(USAGE:
 	$0 -f <filter genelist> -i <input file>
 	[OPTIONAL]
 		-c <config file>
 	
 ) if (!defined $opt_f || ! defined $opt_i);
 my $href;
if ($opt_f){
	$href=getGeneListHash($opt_f);
	if (ref($href)!~/HASH/i){
		print "This failed\n";
	}
}
if (-e $opt_i){
	system ("cp $opt_i $opt_i.orig\n");
	open (INPUT ,"<$opt_i.orig") or die "Cannot read the input file\n";
	open (OUTPUT,">$opt_i") or die "Cannot write to the output file \n";
	my $gene_idx_col=my $impact_idx_col= -1;
	IMPACT: while (my $impact=<INPUT>){
		chomp $impact;
		my @line=split("\t",$impact);
		if ($gene_idx_col<0){
			for (my $i=0;$i<=$#line;$i++){
				if ($line[$i]=~/^Gene$/){
					$gene_idx_col=$i;
				}elsif ($line[$i]=~/^#ANNOVAR annot/i){
					$impact_idx_col=$i;
				}
				last if ($gene_idx_col>-1 && $impact_idx_col>-1);
			}
			print OUTPUT "$impact\n";#hopefully headers
		}else{
			next IMPACT if ($line[$impact_idx_col]=~/intergenic/i);
			my $impact_gene=uc ($line[$gene_idx_col]);
			my @genes=split(/[,:;]/,$impact_gene);
			foreach my $xgene (@genes){
				if ($xgene=~/(\w+)\(/){
					$xgene=$1;
				}
				if (exists $$href{uc($line[$gene_idx_col])}){
					print OUTPUT "$impact\n";next IMPACT;
				}
			}
		}
	}
	close INPUT;close OUTPUT;
}else{
	die "Your file ($opt_i) does not exists in the directory\n";
}