#!/usr/local/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin/";
use ABCC;
use Tie::File;
my @array;
tie @array, 'Tie::File', $ARGV[0] or die "cannot open $ARGV[0] for reading";
if ($ARGV[1]=~/(\d+)\-(\d+)/){
	my $one=$1;my $two=$2;
	print "($one,$two)\t". join("\n",@array[$1..$2])."\n";
}else{
	print $array[$ARGV[1]]."\n";
}
	
