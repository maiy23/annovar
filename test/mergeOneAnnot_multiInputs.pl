#!/usr/bin/perl
use strict;
use Data::Dumper;
#if (`whoami`!~/vuonghm/){die "You are not authorized to use this script\n";}
my $usage= qq(
$0 <annot file> <directory of inputs> 
);
if ($#ARGV!=1){die $usage;}
my $annot_file=$ARGV[0];
my $dir=$ARGV[1];
my $bin=`dirname $0`;chomp $bin;
my $idx_start;my %annot_hash;
open (ANNOT,"<$annot_file" ) || die "Cannot open $annot_file\n";
print "Step 1: Reading in annotation file\n";
while (my $line=<ANNOT>){
	chomp $line;
	my @annot=split("\t",$line);
	if (!$idx_start){
		for (my $i=0;$i<=$#annot;$i++){
			$idx_start=$i if ($annot[$i]=~/^Chr$/);
			last if ($idx_start);
		}
	}
	my $header=join ("\t",@annot[$idx_start..$idx_start+4]);

	$annot_hash{$header}=join ("\t",@annot);
  if ($annot[$idx_start+3] eq '-'){
    $header=join ("\t",@annot[$idx_start..$idx_start+3]);
    $annot_hash{$header}=join ("\t",@annot);
  }
}
close ANNOT;
my @files=`find $dir -maxdepth 2 | grep ANNOVAR.input\$ |sort -u`;
print "Step 2: Finding sample directories\n";
foreach my $dir2run (@files){
#	next if ($dir2run!~/0001/);
  print "Step 3: Working on $dir2run";
  chomp $dir2run;
  my $file_dir=`dirname $dir2run`;chomp $file_dir;
  my $file=`basename $dir2run`;chomp $file;
  open (FILE,"<$dir2run" ) or die "Cannot open $dir2run\n";
  open (OUT,"> $file_dir/$file.annovar_wrpr.output");
  print "\tprinting output to $file_dir/$file.annovar_wrpr.output\n";
  print OUT $annot_hash{"Chr\tQuery Start\tQuery End\tAllele1\tAllele2"}."\n";
  while (my $line=<FILE>){
  	next if ($line=~/^\s+$/);
  	$line=~s/[\r\n]//;
  	my @file_info=split("\t",$line);
  	my $header=join ("\t",@file_info[0..4]);
  	my $info=($#file_info>=5)?join ("\t",@file_info[5..$#file_info]):'';
  	if (exists $annot_hash{$header}){
  		print OUT "$annot_hash{$header}\t$line\n";
  	}elsif ($header=~/#/){
  	}else{
      if ($file_info[3] eq '-'){
        $header=join ("\t",@file_info[0..3]);
        if (exists $annot_hash{$header} ){
          print OUT "$annot_hash{$header}\t$line\n";
          }else{

    		    die "|$header| does not exist($info)$line\n";
          }
        }
  	}
  }
  close FILE;
  close OUT;
  #my $dirname=`basename $dir2run | cut -f1-2 -d"_"`;chomp $dirname;
} 
