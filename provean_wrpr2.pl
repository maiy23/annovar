#!/usr/bin/perl
=head
Written by Hue Vuong
CISB Group
This script is a wrapper to download FASTA files from CyPRUS, generate a variants list and generate a mapper file (gene symbol -> Uniprot Acc) 
and generate bat files to run Provean
When provean has completed, it can also aggregate the results for AVIA protpos annotation (a database file for AVIA)
v1 takes the longest protein as FASTA
#begin mod 2016-06-24
v2 accounts for multiple isoforms by checking the position of the variant and the aa with the fasta sequence
	also organizes subdirectories by alphabet so you can run in batches
=cut
use strict;
use Getopt::Std;
use Data::Dumper;
use lib '/users/abcc/vuonghm/scripts.dir/perl_modules/JSON-2.90/lib/';
use JSON;
use FindBin;
use lib "$FindBin::Bin/";
use ABCC_utils;
use vars qw($opt_i $opt_f $opt_a $opt_d $opt_D $opt_p $opt_m $opt_c $opt_q $opt_e $opt_n);
getopts("f:i:d:D:a:p:c:q:e:n:m");
my $usage=qq(
	$0 
	[REQUIRED]
		-i list of protein coding Variants (ProtPos from AVIA output)
		-a <folder Name> 
			Outputs results into this directory; makes directory if it does not exist
	[OPTIONAL]
		-p <sift, provean> Default : provean
			Certain combinations can be used together
			Acceptable options in [brackets] and the order of operations for each
			[ sift OR provean] - mapping > print mapping file > print fasta > sift/provean > print PBS bat files > execute on PBS
			[ mapping ] mapping >print mapping file
			[ fasta ] mapping > print mapping file > print fasta
			[ executeOnly (sift/provean) ] sift/provean executeOnly
			[ mapback ] mapping > print mapping file > writes AVIA file (\$opt_q) REQUIRED 

			sift | provean | download | mapping | fasta | mapback
			The following must be used in conjuction with (sift|provean|mapback)
			executeOnly
		-D <UCSC build> default mm10;
		-c <filename of previous mapping> If you have run previously and do not need to run here
			INPUT -> OUTPUT
		-e <filename of previous mapping that has data that you wish to EXCLUDE> If you have run previously and do not need to run here
		-q <ORIGINAL file (not just exonic) -> This option is specific for AVIA
		-n <output name for mapback only>
);
if ($opt_p !~/mapback/ && (!defined $opt_i || !defined ($opt_a))){
	die "You must specify an input file ($opt_i) and a results directory($opt_a)\n\n$usage\n";
}
if ($opt_a!~/^\//){
	die "Please specify the full path the output directory\n";
}
my $debug=1;
$opt_p=(defined $opt_p)?lc($opt_p):"provean";
$opt_D||="mm10";
my $taxon=getOrgType("$opt_D","taxonid");
##Paths to scripts and executables

## for AVIA users, we expect the input ($opt_i) to the converter script to look like:
# NM_001011874:D610D,
# NM_001011874:P528P,
# NM_001011874:G514G,
# NM_001011874:F510F,
# NM_001011874:I408I,
# NM_001011874:I371I,
my $converterScript="https://bioinfo-abcc.ncifcrf.gov/cyprus/rest/v1/getSequences?taxon=$taxon&ids";#needs a file added
my $mergeScript='';#this merges the output of one into the other
my $sift4G_exe='';
my $provean_exe='/opt/nasapps/development/provean/1.1.5/bin/provean.sh ';
my $suffix;
if ($opt_p=~/provean/i){
	$suffix='.var';
}elsif ($opt_p=~/sift/i){
	$suffix="-X.subst";
}else{
	$suffix=".var";
}
if ($debug){
	open (LOG,"|tee $opt_a/log") or die "cannot open $opt_a/log";
}else{
	open (LOG,">$opt_a/log") or die "Cannot open $opt_a/log for writing\n";
}
## Set up
mkdir("$opt_a") if (!-e $opt_a);
my $delimiter='\t';
#Map to the input file
#putting FASTA into separate files for each protein (take longest UniProt Acc)
#and the vars in another file in the specified users directory 
my %mapping;my @list;
if (!defined $opt_c || !-e $opt_c){
	##previous mapping file gene->uniprotID not found, run this block script
	die "Cannot find ($opt_i) for reading at line ". __LINE__."\n" if (!-e "$opt_i");
	print LOG "[INFO] Reading your input $opt_i\n" if ($debug);
	open (INPUT,"<$opt_i") or die "Cannot open $opt_i for reading\n";
	while (<INPUT>){
		my @trxArr=split(",",$_);chomp $trxArr[$#trxArr];
		foreach my $trxLine (@trxArr){
			my ($id,$protPos)=split(":",$trxLine);
			chomp $protPos;
			if ($protPos=~/(p.){0,1}([A-Z]\d+[A-Z])/){
				$protPos=$2;
			}else{
				print LOG "skipping ||$protPos||\n";
				next;
			}
			push(@{$mapping{$id}},$protPos);
		}
	}
}elsif ($opt_p=~/mapback/){
	#Read in the input file to hash\
	open (MAPPER,"<$opt_c") or die "Cannot read mapping file $opt_c";
	my @array=<MAPPER>;
	foreach my $line (@array){
		my ($id,$mapped)=split("\t",$line);chomp $mapped;
		$mapping{$mapped}=$id;
	}
}else{
	print LOG $usage and exit;
}
my $printMapper='';
my $count=0;
print LOG "[INFO] About to run mapping using CyPRUS\n" if ($debug);

my $genecount=0;my $batchGenes='';my $res;
##First run CyPRUS in batches
if ($opt_p=~/(fasta|sift|provean|mapping)/){
	foreach my $gene (sort keys %mapping){
		next if (!$gene);
		if ($genecount>=50){
			# print "Change this back to 50 at line" . __LINE__ . "\n";
			run();
			$batchGenes='';$genecount=0;
		}
		$batchGenes.="$gene,";$genecount++;
	}
	if ($genecount){
		print LOG "Running last batch : $genecount\n";
		run();
	}
	print LOG "[INFO] Done! mapping using CyPRUS\n" if ($debug);
}


##Now print mapback file
if (!defined $opt_c){#This is used to map back 
	print LOG "[ERR] Could not write to mapper file ($opt_i.mapper) at line ". __LINE__."\n" and exit if (printTo(">$opt_i.mapper", $printMapper));
}

exit if ($opt_p=~/(fasta|mapping)/);
print LOG "[INFO] About to $opt_p...\n" if ($debug);
if ($opt_p=~/provean/i){
	opendir(my $dh, "$opt_a") or die "Cannot open dir $opt_a\n";
	while (readdir $dh){
		if ($_=~/(.*).fasta$/){
			my $base=$1;
			print printPBS ("$provean_exe -q $opt_a/$base.fasta -v $opt_a/$base$suffix > $opt_a/$base.results","$opt_p\_$base","$opt_a") ."\n";
		}
	}
}elsif ($opt_p=~/sift/i){
	#sift4G can have all files in a single directory
	#run the entire directory
	die 'Havent fully coded this part but you could probably add this in here at line ' . __LINE__. "\n";
}

## Finally, we aggregate the files for AVIA
if($opt_p=~/mapback/){
	open (DIR,"$opt_a") or die "Cannot open $opt_a for directory walking\n";
	$opt_n||="$opt_D\_$opt_p.txt";
	open (OUT,">$opt_n") or die "Cannot write to $opt_n";
	my @dirs=split("\n",`ls $opt_a`);
	foreach my $subdir (@dirs){
		chomp $subdir;
		if (-d "$opt_a/$subdir"){#walk the directory for original run since it was very large
			open (SUBDIR,"<$opt_a/$subdir") or die "cannot open $opt_a/$subdir\n";
			my @subdirs=split("\n",`ls $opt_a/$subdir`);
			foreach my $resfn (@subdirs){
				chomp $resfn;
				mapback("$opt_a/$subdir/$resfn") if ($resfn=~/results$/);
			}
			close SUBDIR;
		}elsif ($subdir=~/.results$/){#if the files are already in the directory, then evaluate them
			mapback("$opt_a/$subdir",\*OUT);	
		}
	} 
	close DIR;close OUT;
}

sub mapback{
	my $fn=shift;
	open (FILE,"<$fn") or die "Cannot open $fn\n";
	my $uni_id=`basename $fn`;chomp $uni_id;
	$uni_id=~s/.results//;
	my $avia_gene; 
	if (exists $mapping{$uni_id}){
		$avia_gene=$mapping{$uni_id};
	}else{
		warn "$uni_id does not exist in mapping!\n";
		return;
	}
	while (<FILE>){
		if ($_=~/^([A-Z]\d+[A-Z])\t(\S+)/){
			my $var=$1;my $score=$2;
			print "$avia_gene:$var\t$score\n";
		}
	}
	close FILE;
}


sub run{
	#now do the mapping\
	print LOG "ABout to run curl -s '$converterScript=$batchGenes'\n";
	$batchGenes=~s/,$//;
	$res=(`curl -s '$converterScript=$batchGenes' `);
	if ($res){
		$res=decode_json($res);
	}else{
		print LOG "Skipping $batchGenes b/c CyPRUS couldn't find it(empty)\n";return 0;
	}
	if (ref($res) eq 'ARRAY' ){ #CyPRUS returns an array when empty
		print LOG "Skipping $batchGenes b/c CyPRUS couldn't find it\n";next GENE;
	}
	GENE: foreach my $gene (sort keys %{$res}){
		next if (!$gene);
		chomp $gene;
		my $longest=0;my $longestId='';my $uniprot_ct=0;
		foreach my $unId (keys $$res{$gene}){
			my $length=$$res{$gene}{$unId}{'length'};
			$longestId=$unId and $longest=$length if ($length > $longest);
			my $sequence=$$res{$gene}{$unId}{'sequence'};
			##validate that all variants will match reference
			my @matchingVars;
			VAR : foreach my $var (@{$mapping{$gene}}){
				my $pos = my $aa;
				if ($var=~/(\w)(\d+)(\w)/){
					$pos=$2;$aa=uc($1);my $other=$3;
					# print "skipping $var\n" and next VAR if ($aa eq $other);
				}
				if ($pos> $length){
					next VAR;
				}
				my $prot_aa=uc(substr($sequence,($pos-1),1));
				if ($prot_aa eq "$aa"){#check that it matches
					push(@matchingVars,$var);
					print "\t$prot_aa is EQUAL to $aa!\n";
					# $found{$var}=1;
				}
			}
			if ($#matchingVars>-1){
				if ($opt_p=~/(fasta|provean|sift)/ && $opt_p!~/executeOnly/){
					my $fasta=formatFasta($sequence);
					print LOG "About to write to $unId.fasta in $opt_a\n" if ($debug);
					open (FASTA,">$opt_a/$unId.fasta") or die "Cannot open ($unId.fasta)\n";
					print  FASTA ">$unId\n$fasta";
					close FASTA;
					open (SNPS,">$opt_a/$unId$suffix") or die "Cannot open ($unId$suffix)\n";
					print SNPS join ("\n",@matchingVars);
					close SNPS;
				}
				$printMapper.="$gene\t$unId\n";
			}
			#delete this for subsequent runs
			# we did this because we used the longest ones originally
			# now we are trying to catch all the variants that we missed
			if (-e "$opt_a/$longestId.fasta"){
				print "deleting $longestId\n";
				unlink ("$opt_a/$longestId.fasta");
				unlink ("$opt_a/$longestId.var");
			}
		}
	}
	return 1;
}