#!/usr/local/bin/perl 
#use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
umask(0000);
use Data::Dumper;
=head
This module describes ABCC/LMT specific modules for additional calculations not done by ANNOVAR
v2.0
bug fixes to the output where errors occur
=cut
use FindBin;
use lib "$FindBin::Bin/";
use GSMutate;

my $fa=getmRNASeq('NM_002749','');
print $fa;