#!/usr/local/bin/perl 
#use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
umask(0000);
my $dbPath='/SeqIdx/annovardb/';
use Data::Dumper;
=head
This module describes ABCC/LMT specific modules for additional calculations not done by ANNOVAR
v2.0
bug fixes to the output where errors occur
=cut

sub getmRNASeq{
	my $id=shift;
	my $seq_fn=shift;
	#Set up defaults
	if (!$seq_fn || $seq_fn=~/refgene/){
		$seq_fn=$dbPath."humandb/hg19_refGeneMrna.fa";
	}elsif($seq_fn=/ensembl/){
		$seq_fn=$dbPath."humandb/hg19_ensGeneMrna.fa";
	}

	if(!-e $seq_fn){
		die "Your file does not exist ($seq_fn)";
	}

	$fa=`grep '>$id' -A 1 $seq_fn`;
	return $fa;
}
sub getGeneSequence{
	my $geneName=shift;
}
sub mutate{
	my $seq=shift;
	my $mt=shift;
	my @arr=split("",$seq);



}
sub runGenScan{
	
}
1;