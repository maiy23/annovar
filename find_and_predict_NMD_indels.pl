#!/usr/bin/perl
use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use vars qw( $opt_o $opt_c $opt_d $opt_i $opt_I $opt_v);
getopts("c:d:i:o:I:v:");
#GLOBAL DEPENDENT VARS
my $SIFT_exome_cmd='/bioinfoA/sift/sift_mysql/bin/SIFT_exome_indels.pl';
my $CDS_DIR='/SeqIdx/siftdb/coding_info_37';
my $VARIATIONS_DIR='/SeqIdx/siftdb/Human_db_37';

if ($opt_c && -e $opt_c){
	$CDS_DIR=$opt_c;
}
if ($opt_d){
	$VARIATIONS_DIR=$opt_d;
}
my $usage=qq(This script was written for the AVIA annovar pipeline and requires a very specific input format
	You must specify an ANNOVAR output file!
	$0
	[REQUIRED] 
		 -i <ANNOVAR.wrpr.output file>
		 -I <not collapsed
	[OPTIONAL]  #defaults are set internally
		 -c <full path to the SIFT's coding info files>
		 -d <name of the SIFT databasename> 
		 -o <output name> DEFAULT yourinputname.IndelPred
);
die "$usage\n" if (!defined $opt_i && !defined $opt_I);
my $input=(defined $opt_i)?$opt_i:$opt_I;
$opt_o||="$input.IndelPred";
die "Your file is empty or does not exist($input)\n" if (! -e "$input" || -z $input);
my $dir='./';
if ($input=~/^\//){
	$dir=`dirname $opt_i`;my $basename=`basename $opt_i`;chomp $dir;chomp $basename;
	chdir($dir);
	$opt_i=$basename;
}

system ("grep -e fs -e indel -e ins -e del $input -n | grep exonic | grep -ve ncrna -i  >$input.exonic_indels.tmp\n");
if (-z "$input.exonic_indels.tmp"){system ("rm $input.exonic_indels.tmp\n");exit;}
open (TMP,"<$input.exonic_indels.tmp") or die "Cannot open file\n";
open (FINAL,"> $input.exonic_indels.csv") or die "cannot write to file\n";
my $idx=-1;my %indel;
if ($opt_i){
	while (my $line=<TMP>){
		chomp $line;
		my @arr=split("\t",$line);
		if ($idx==-1){
			until ($arr[$idx]=~/exonic/ || $idx>$#arr){
				$idx++;
			}
			die "could not find exonic index" if ($idx>=$#arr);
		}
		#these positions are all relative to the ANVR annot flag all standard AVIA output workflows
		if ($arr[7+$idx] eq '-' || $arr[8+$idx] eq '-'){#check that it is truly an indel
			my $line=($arr[0]=~/^(\d+):/)?$1:'';
			$arr[4+$idx]=~s/chr//g;
			if ($arr[8+$idx] eq '-'){#deletion
				print FINAL join (",",$arr[4+$idx],$arr[5+$idx],$arr[6+$idx]+1,$arr[3+$idx]."1","$arr[7+$idx]/$arr[8+$idx]",$line)."\n";
				
			}else{
				print FINAL join (",",$arr[4+$idx],$arr[5+$idx],$arr[6+$idx],$arr[3+$idx]."1","$arr[7+$idx]/$arr[8+$idx]",$line)."\n";
			}
			$indel{$line}='';
		}
	}
}elsif ($opt_I){
	while (my $line=<TMP>){
		chomp $line;
		my @arr=split("\t",$line);
		if ($idx==-1){
			until ($arr[$idx]=~/exonic/ || $idx>$#arr){
				$idx++;
			}
			die "could not find exonic index" if ($idx>=$#arr);
		}
		#these positions are all relative to the ANVR annot flag all standard AVIA output workflows
		if ($arr[7+$idx] eq '-' || $arr[8+$idx] eq '-'){#check that it is truly an indel
			my $line=($arr[0]=~/^(\d+):/)?$1:'';
			$arr[4+$idx]=~s/chr//g;
			if ($arr[8+$idx] eq '-'){#deletion
				print FINAL join (",",$arr[4+$idx],$arr[5+$idx],$arr[6+$idx]+1,$arr[3+$idx]."1","$arr[7+$idx]/$arr[8+$idx]",$line)."\n";
				
			}else{
				print FINAL join (",",$arr[4+$idx],$arr[5+$idx],$arr[6+$idx],$arr[3+$idx]."1","$arr[7+$idx]/$arr[8+$idx]",$line)."\n";
			}
			$indel{$line}='';
		}else{
			print "not an indel $arr[7+$idx] // $arr[8+$idx]\n";
		}
	}
}
close TMP;close FINAL;
#run the sift exome script
## use this perl version for moab b/c of SQLite.pm and Dynaloader
# die "/opt/nasapps/applibs/perl-5.16.2/bin/perl $SIFT_exome_cmd -i $input.exonic_indels.csv -c $CDS_DIR -o $dir -d $VARIATIONS_DIR 2>/dev/null\n";
print STDERR "[INFO] Running script: /opt/nasapps/applibs/perl-5.16.2/bin/perl $SIFT_exome_cmd -i $input.exonic_indels.csv -c $CDS_DIR -o $dir -d $VARIATIONS_DIR 2>/dev/null\n";
my $pid=`/opt/nasapps/applibs/perl-5.16.2/bin/perl $SIFT_exome_cmd -i $input.exonic_indels.csv -c $CDS_DIR -o $dir -d $VARIATIONS_DIR 2>/dev/null`;chomp $pid;
# my $pid='Your job id is 28710 and is currently running.';
if ($pid=~/Your job id is (\d+) and is currently running./){
	$pid=$1;
}else{
	die "could not parse this ($pid)\n";
}
open (PRED,"<$pid/$pid\_predictions.tsv") or die "cannot open the output files or it does not exist $pid/$pid\_predictions.tsv\n";
while (my $line=<PRED>){
	chomp $line;
	#1based (not 0)column 10 is the Causes NMD, col 12 is the transcript viz,col 7 is the % indel location
	 my @line_arr=split("\t",$line);
	 if ($line_arr[0]=~/,(\d+)\s*$/){
	 	my $pred;my $lineid=$1;
	 	if ($line_arr[9]=~/yes/i){
	 		$pred="Causes NMD";
	 	}else{
	 		$pred="Does not cause NMD";
	 	}	 	
	 	$indel{$lineid}=join (":",$pred,$line_arr[6]);#pos 11 is the nice picture
	 }else{
	 	die "$line_arr[0]|" if ($line!~/Coordinate/);
	 }
}
close PRED;

open (NEW,">$opt_o") or die "Cannot open output for writing\n";
$idx='NOTFOUND';#means not found,reused
my $linect=0;
if ($opt_i){
	open (FILE,"<$input") or die "Cannot open file $input for reading";
}elsif (defined $opt_I){
	my $newOut=$input;
	$newOut=~s/.variant_function.collapsed/.$opt_v/;
	open (FILE,"<$newOut") or die "Could not find $0:".__LINE__.":$newOut for reading\n";
	print "reading $newOut\n";
	$idx=0;
}
while (my $sift_line=<FILE>){
	$linect++;
	chomp $sift_line;
	if ($sift_line=~/^#/ || $idx eq 'NOTFOUND'){
		#find the predictions
		my @line=split("\t",$sift_line);
		for (my $i=0;$i<=$#line;$i++){
			if ($line[$i]=~/sift/i){
				$idx=$i;
			}
			last if ($idx=~/\d+/);
		}
		die "Could not find the index of the sift scores\n" if ($idx eq 'NOTFOUND');
	}
	if (exists $indel{$linect}){
		my @line=split("\t",$sift_line);
		$line[$idx]=$indel{$linect};
		print NEW join ("\t",@line)."\n";
	}else{
		print NEW "$sift_line\n";
	}
}
close FILE;
close NEW;
#delete intermediate files
system ("rm $input.exonic_indels.tmp $input.exonic_indels.csv \n");
