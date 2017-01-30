#!/usr/bin/perl
use strict;
# Read $ARGV[1] mapping file (Gene Syn -> UniProt Acc) and 
# walk through a directory with $ARGV[0] suffix Default .results
# mapping GeneSym to UniProt and printing out file

##DEFAULTS;
my $suffix='.results';
my $mapping_fn='mouse.coding.ProtPos.mapper';
my $directory='./';
my $output="mm10_provean.txt";

if ($#ARGV>-1){
	$mapping_fn=$ARGV[0];
}
if($#ARGV>0){
	$suffix=$ARGV[1];
}


print "[INFO] Mapping $mapping_fn to all results in *$suffix in the current directory\n";
my %mapper;
open (MAPPER,"<$mapping_fn") or die "Cannot open $mapping_fn\n";
while (<MAPPER>){
	my ($genesym,$uniprot_acc)=split("\t",$_);chomp $uniprot_acc;
	$mapper{$uniprot_acc}=$genesym;
}
close MAPPER;
open (OUT,">>$output") or die "Cannot open $output for writing\n";
# now walk through the directory and write to an output file
opendir (DIR, "./") or die "Cannot open ./\n";
my @files = readdir(DIR);
closedir(DIR);
foreach my $file (@files) {
	if (($file !~ /^\.\.?$/) && $file=~/$suffix$/) {
		my $name;
		if ($file=~/(\S+).results/){
			$name=$1;
		}else{
			die "Could not parse $file (\\S+).results\n";
		}
		my $gene;
		if (exists $mapper{$name}){
			$gene=$mapper{$name};
		}else{
			die "$name does not exist in the mapping file\n";
		}
		print "reading $file\n";
		open (FH,"<$file") or die "cannot open $file\n";
		while (<FH>){
			next if ($_=~/^(#|\[)/);
			my ($variation,$score)=split("\t",$_);chomp $score;
			if ($score<=-2.5){
				$score="Deleterious:$score";
			}else{
				$score="Neutral:$score";
			}
			print OUT "$gene:$variation\t$score\n";
		}
		close FH;
	}
}
print "Done!\n";
