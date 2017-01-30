#!/usr/bin/perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/current";
use ABCC;
use Getopt::Std;
use File::stat;
use Time::localtime;
use vars qw( $opt_f $opt_o $opt_d $opt_u $opt_i);
 getopts("d:o:f:ui");
umask(0000);
#################### Set up global directories and scripts ##############

#by default the directory is /bioinfoC/AVA/INTERNAL_SUBMISSIONS
$opt_d||='/bioinfoC/AVA/INTERNAL_SUBMISSIONS';
$opt_o||='/bioinfoC/AVA/PROC';
our $bin='/bioinfoC/hue/annovar/current';
#########################################################################
#read the directory
#match up sets of data
my @files=`find $opt_d -name '*.txt' | grep -e Blood -e Tumor`;
my %paired;
my $curr_time=time();
foreach my $fullpath (@files){
	chomp $fullpath;
	my $mfile = (stat $fullpath);
	my @structure=split("/",$fullpath);
	my $file=pop@structure;
	my $diff=$curr_time-$$mfile[9];
	my $sample;
	if ($file=~/(OV\_\d+)\_[a-z]+(\_[a-z]){0,1}\_Def(Snps|Dips)\.txt/i){
		$sample=$1;
		push (@{$paired{$sample}{files}},$fullpath);
	}
	if ($opt_i || ! exists $paired{$sample}{time} || (exists $paired{$sample}{time} && $paired{$sample}{time}<$diff)){
		$paired{$sample}{time}=$diff;
	}
}
my $runFilter=0;my %files;
my $count=0;
foreach my $sample (keys %paired){
#	$sample='OV_0001';
	next if (!$opt_u && $#{$paired{$sample}{files}}<3);#waiting for 4 files
	my $filedir=($#{$paired{$sample}{files}}<3)?"$sample\_partial":$sample;
	%files=();
	next if (!$sample);
	print "$sample is too young\n" and next if ($paired{$sample}{time}< 1200);#too young 20 minutes because it could be still copying
	next if (-e "$opt_o/$filedir");
	$count++;
	#make a sample specific directory under PROC
	mkdir ("$opt_o/$filedir");
	chdir ("$opt_o/$filedir");
	SAFE_PROCESS ("cp $opt_o/config.ini .\n");
	foreach my $file (@{$paired{$sample}{files}}){
		chomp $file;my $aa;#=`basename $file`;chomp $aa;
		if ($file=~/\_BL\_/){
			$aa="hg19_$sample\_Blood.txt";
		}else{
			$aa="hg19_$sample\_Tumor.txt";
		}
		#convert all 
		SAFE_PROCESS ("perl $bin/convert2annovar.pl --format casava --af --outfile $aa --append --header --logfile $aa.cvrt2anvr.log $file \n");
		$files{$aa}=1;
	}
	if ($runFilter){
		open (DB,">searchTheseDbs.txt") or die "Cannot open searchTheseDbs.txt\n";
		print DB "$opt_o/$filedir/hg19_$sample\_Blood.txt\n";
		print DB "$opt_o/$filedir/hg19_$sample\_Tumor.txt\n";
		close DB;
		SAFE_PROCESS("cat $opt_o/searchTheseDBs.txt >> searchTheseDbs.txt\n");
		my $success=_index("$opt_o/$filedir","hg19_ANVR_Blood.txt");
		if (!$success){
			die "Could not index your file!\n";
		}
		#do subtractive analysis A-B and B-A
		SAFE_PROCESS("perl $bin/annotate_variation_ABCC.pl --filter --dbtype generic --genericdbfile hg19_$sample\_Blood.txt hg19_$sample\_Tumor.txt .\n");
		SAFE_PROCESS("perl $bin/annotate_variation_ABCC.pl --filter --dbtype generic --genericdbfile hg19_$sample\_Tumor.txt hg19_$sample\_Blood.txt .\n");
		SAFE_PROCESS ("ln -s hg19_$sample\_Tumor.txt.hg19_$sample\_Blood.exact_filtered ANNOVAR.input\n");#Tumor minus normal
		SAFE_PROCESS ("ln -s hg19_$sample\_Blood.txt.hg19_$sample\_Tumor.exact_dropped LOH.input\n");#normal minus tumor
		#set up the .abcc web file or config file
		SAFE_PROCESS ("cat $opt_o/searchTheseDBs.txt >> searchTheseDbs.txt;cp $opt_o/config.ini .\n");
		#run processing_wrapper
		SAFE_PROCESS( "perl $bin/annovar_qsub_wrpr.pl -f searchTheseDbs.txt -x -z -R -i ANNOVAR.input\n");
		
	}else{
		foreach my $aa (keys %files){
			next if ($aa=~/normal/i);
			open (DB,">$aa.searchTheseDbs.txt") or die "Cannot open searchTheseDbs.txt\n";
			print DB "$opt_o/$sample/$aa.hg19_af\n";
			close DB;
			SAFE_PROCESS("cat $opt_o/searchTheseDBs.txt >> $opt_o/$filedir/$aa.searchTheseDbs.txt\n");
			SAFE_PROCESS( "perl $bin/annovar_qsub_wrpr.pl -f $aa.searchTheseDbs.txt -x -z -R -i $aa \n");
		}
	}
	exit if ($count==4);
}
	
__DATA__
user.label=FOO
user.file=FOO
user.internal=1
user.modulesid=3
user.date=FOOf
user.basedir=/bioinfoC/AVA/PROC
ref.ver=hg19
user.inputformat=casava
user.calcAF=true
user.modulesid=3
user.name_column_header=1
user.show_flanking=1
report.annotdb_siftv63=on
report.annotdb_ljb_pp2=on
report.annotdb_cosmicdb=on
report.annotdb_snp135=on
report.annotdb_cg69=on
report.annotdb_hapmap_3_3=on
report.annotdb_1000gALL_sites_2012_04=on
report.annotdb_gwasCatalog=on
report.annotdb_targetScanS=on
report.annotdb_tfbsConsSites=on
user.program=viz
