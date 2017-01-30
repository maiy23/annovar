#!/usr/bin/perl
use strict;
use Data::Dumper;
umask(0000);
use Getopt::Std;
use vars qw( $opt_i $opt_o $opt_k $opt_E );
getopts("i:o:k:E");#database and avia id, respectively
my $input=$opt_i;
my $outs=(defined $opt_o )?$opt_o:"$input.out";
my $donotcheck=(defined $opt_k)?$opt_k:1;
my $useensembl=(defined $opt_E)?"Ensembl_":'';
die "You must supply a valid annotation file\n" if (!$input || !-e $input);
open (FILE,"<$input") || die "Cannot open $input\n";
open (OUT,">>$outs") || die "Cannot open $outs\n";
my ($end,@headers);
my %hash;#unique hash;
while (my $line=<FILE>){
  chomp $line;
  my @arr=split("\t",$line);chomp $arr[$#arr];
  if ($line=~/^#/){
	$line=~s/#//g;
	@headers=split("\t",$line);
	## find the column in which the ANNOVAR input was found
	if ($line=~/(^.*\tChr)\tQuery Start\tQuery End\tAllele1\tAllele2/){
		my $string=$1;
		$end=0;
		while ($string =~ /\t/g) { $end++ }
	}else{
		die "could not find (^.*\tChr)\tQuery Start\tQuery End\tAllele1\tAllele2 in ($line)\n";
	}
  }else{
	die "Could not find your header row\n" if (!$end || $end<0);
	if (length($arr[$end+3]) > 100){$arr[$end+3]=".";}#change the reference allele to "." for large indels
	my $pos=join ("\t",@arr[$end..$end+4]);
	die join ("\n",@headers). "\n\tEND=$end\n\t$pos" if ($pos!~/^(chr){0,1}[\dXYMT]{1,}\t\d+\t\d+/);
	for (my $i=0;$i<=$end-2;$i++){
		my $db_row=$arr[$i];

  		if (!$donotcheck && exists $hash{"$pos\t$headers[$i]"}){next;}
		if (length($arr[$i])>2000){
			my $substring = substr ($arr[$i],0,1000)."...";
			#print "substringing($arr[$i]) $substring\n";
			$arr[$i]=$substring;
		}
		if ($arr[$i] eq 'ERR' || $arr[$i]=~/INTERNAL ERR/i || $arr[$i] eq ''){

		}elsif ($i>$end -5){
			print OUT "$pos\t$useensembl$headers[$i]\t$arr[$i]\n";
			$hash{"$pos\t$headers[$i]"} if (!$donotcheck);
		}else{
			print OUT "$pos\t$headers[$i]\t$arr[$i]\n";
			$hash{"$pos\t$headers[$i]"} if (!$donotcheck);
		}#write empty if '-'
	}
  }
}
close OUT;close FILE;
system ("chmod 777 $outs\n");
