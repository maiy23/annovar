#!/usr/local/bin/perl
use strict;
use Data::Dumper;
umask(0000);
my $infofile='/SeqIdx/annovardb/humandb/database.info';
if (!-e $infofile){die "$infofile no longer exists!! please generate\n";}
die "Please enter a filename\n " if ($#ARGV==-1);
open (FILE,"<$ARGV[0]") or die "Cannot open file for reading";
my %headers;my @final_headers;
open (DB,">avia.database.info.txt") or die "Cannot open avia.database.info.txt for writing\n";
my $change_header=0;
if ($#ARGV==0){
	$change_header=1;
	open (OUT,">$ARGV[0].rename.tmp") or die "cannot write $ARGV[0].rename.tmp for writing\n";
}else{
	open (OUT,">subtractive_headers.out") or die "cannot open subtractive...for writing\n";
}
my $headerfound=0;
while (<FILE>){
	if ($_=~/^(#|Chr)/i && $_=~/(#{0,1}ANNOVAR\sannot|ANNOVAR_ANNOT)/i){
		$headerfound++;
		my @header=split("\t",$_);chomp $header[$#header];
		if ($change_header){
			print DB "#Database\tDatabase Description\tSource\tVersion\n";
		}else{
			print DB join ("\t",'AVIA','Database','Database Description','Source','Version')."\n";
		}
		for (my $i=0;$i<=$#header;$i++){
			$header[$i]=~s/^(#|\s)//g;
			print OUT "#$header[$i]\t" and next if ($header[$i]=~/^(\s*|Summary)$/);
			# last if ($header[$i]=~/(ANNOVAR annot)/i && !$change_header);
			print OUT "#$header[$i]\t" and next if ($header[$i]=~/(FlankingSequence|samplename)/i && !$change_header);
			$headers{$header[$i]}=$i;
			# if ($change_header){print OUT join ("\t",@header[$i..$#header])."\n" and last if ($header[$i]=~/ANNOVAR\sannot/i);}
			print OUT "#$header[$i]\t" and next if ($header[$i]=~/(ANNOVAR(\s|\_)annot|Variant ID|^\s*Gene$|(ANNOVAR|Annot)(\s|\_)Feat|Rel pos|Gene(\s|\_)Ori|Query.*|Allele\d+|Comment|^Chr$|ProtPos)/i);
			$header[$i]=~s/\.exact$//;
			my @line=split("\t",`grep '\\b$header[$i]\\b' $infofile -m 1`);

			chomp $line[$#line] if ($#line>-1);	
			if ($#line==-1){
				print OUT "#$header[$i]\t" ;
				print DB "$header[$i]\t-\t-\t-\n" ;
			}else{
				chomp $line[$#line];
				print OUT "#$line[2]\t" ;
				my $addon="$line[1]\t" if (!$change_header);
				print DB "$addon$line[2]\t$line[6]\t$line[3]\t$line[4]\n";
			}
			
		}
		print OUT "\n" and print "done!" and close DB and exit if ($#ARGV>0);
		print OUT "\n";
	}else{
		die "Could not find your header!!\n" if ($headerfound==0);
		print OUT $_;
	}
}
if (keys %headers == 0){
	die "Could not find your header\n";
	system ("rm $ARGV[0].rename.tmp\n");
}
close FILE;close DB;close OUT;
if ($change_header){
	system ("mv $ARGV[0].rename.tmp $ARGV[0]\n");
}


