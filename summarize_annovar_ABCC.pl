#!/usr/bin/perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
## ABCC GLOBAL VARIABLES
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use ABCC;#ABCC
our ($keepline,$relpos,$header,$allupdated,$alt_input);#filter vs annotate
### end ABCC VARs
our $VERSION = 			'$Revision: 488 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-10-02 21:27:43 -0700 (Sun, 02 Oct 2011) $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $buildver, $step, $checkfile, $remove, $verdbsnp);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'buildver=s'=>\$buildver, 'step=s'=>\$step, 
	'checkfile!'=>\$checkfile, 'remove'=>\$remove, 'verdbsnp=i'=>\$verdbsnp , 'keepline'=>\$keepline, 'relpos'=>\$relpos, 'header'=>\$header , 'allupdated'=>\$allupdated) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

($queryfile, $dbloc) = @ARGV;

$outfile ||= $queryfile;
$buildver ||= 'hg18';
$buildver eq 'hg18' or $buildver eq 'hg19' or pod2usage ("Error in argument: the --buildver argument can be 'hg18' and 'hg19' only");
not defined $checkfile and $checkfile = 1;

defined $verdbsnp or $verdbsnp = 130;
##abcc
my $abcc_addons='';
if ($header){
	$abcc_addons.=" --header ";
}
if ($relpos){
	$abcc_addons.=" --relpos ";
}
if ($keepline){
	$alt_input= "$queryfile.exon.input";
	$abcc_addons.=" --keepline ";
}
##end abcc
if ($allupdated){#this is where you keep all the up-to-date versions of databases
	$buildver='hg19';
	$verdbsnp='132';
}
my %valistep;
if ($step) {
	my @step = split (/,/, $step);
	for my $i (0 .. @step-1) {
	if ($step[$i] =~ m/^(\d+)-(\d+)$/) {
		for my $nextstep ($1 .. $2) {
		$valistep{$nextstep}++;
		}
	} elsif ($step[$i] =~ m/^(\d+)$/) {
		$valistep{$1}++;
	} else {
		pod2usage ("Error: invalid -step argument ($step) is specified. Please use comma-separated number only (dash line such as 1-5 is accepted)");
	}
	}
} else {
	for my $nextstep (1..12) {
	$valistep{$nextstep}++;
	}
}

$checkfile; # and checkFileExistence ();


my $sc;

#run step 1
if ($valistep{1}) {
	$sc = "annotate_variation_ABCC.pl -geneanno -buildver $buildver -outfile $outfile $abcc_addons -exonsort $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 1 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step2
if ($valistep{2}) {
	if ($buildver eq 'hg18') {
		$sc = "annotate_variation_ABCC.pl -regionanno -dbtype mce44way -buildver $buildver $abcc_addons -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 2 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	} elsif ($buildver eq 'hg19') {
		$sc = "annotate_variation_ABCC.pl -regionanno -dbtype mce46way -buildver $buildver $abcc_addons -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
}

#run step3
if ($valistep{3}) {
	$sc = "annotate_variation_ABCC.pl -regionanno -dbtype segdup -buildver $buildver $abcc_addons -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 3 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}


if ($buildver eq 'hg19') {
	if ($valistep{4} or $valistep{5} or $valistep{6}) {
		$sc = "annotate_variation_ABCC.pl -filter -dbtype ALL.sites.2010_11 -buildver $buildver $abcc_addons -outfile $outfile $queryfile $dbloc";
		if ($keepline){
			$sc=~s/\-filter/\-regionanno/;
		}
		print STDERR "\nNOTICE: Running step 4/5/6 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
} 
else {
	#run step4
	if ($valistep{4}) {
		$sc = "annotate_variation_ABCC.pl -filter -dbtype 1000g2010jul_ceu -buildver $buildver $abcc_addons -outfile $outfile $queryfile $dbloc";
		if ($keepline){
			$sc=~s/\-filter/\-regionanno/;
		}
		print STDERR "\nNOTICE: Running step 4 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
	#run step5
	if ($valistep{5}) {
		$sc = "annotate_variation_ABCC.pl -filter -dbtype 1000g2010jul_yri -buildver $buildver $abcc_addons -outfile $outfile $queryfile $dbloc";
		if ($keepline){
			$sc=~s/\-filter/\-regionanno/;
		}
		print STDERR "\nNOTICE: Running step 5 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
	#run step6
	if ($valistep{6}) {
		$sc = "annotate_variation_ABCC.pl -filter -dbtype 1000g2010jul_jptchb -buildver $buildver $abcc_addons -outfile $outfile $queryfile $dbloc";
		if ($keepline){
			$sc=~s/\-filter/\-regionanno/;
		}
		print STDERR "\nNOTICE: Running step 6 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
}

#run step7
if ($valistep{7}) {
	$sc = "annotate_variation_ABCC.pl -filter -dbtype snp$verdbsnp -buildver $buildver $abcc_addons -outfile $outfile $queryfile $dbloc";
	if ($keepline){
		$sc=~s/\-filter/\-regionanno/;
	}
	print STDERR "\nNOTICE: Running step 7 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step8
if ($valistep{8}) {
	$sc = "annotate_variation_ABCC.pl -filter -dbtype avsift -buildver $buildver -sift 0 $abcc_addons -outfile $outfile $alt_input $dbloc";
	if ($keepline){
		$sc=~s/\-filter/\-regionanno/;
		$sc=~s/\-sift 0//g;
	}
	print STDERR "\nNOTICE: Running step 8 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step9
if ($valistep{9}) {
	$sc = "annotate_variation_ABCC.pl -filter -dbtype ljb_pp2 -score_threshold 0 -buildver $buildver -outfile $outfile $alt_input $dbloc";
	if ($keepline){
		$sc=~s/\-filter/\-regionanno/;
	}
	print STDERR "\nNOTICE: Running step 9 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step10
if ($valistep{10}) {
	$sc = "annotate_variation_ABCC.pl -filter -dbtype ljb_phylop -score_threshold 0 -buildver $buildver $abcc_addons -outfile $outfile $alt_input  $dbloc";
	if ($keepline){
		$sc=~s/\-filter/\-regionanno/;
	}
	print STDERR "\nNOTICE: Running step 10 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step11
if ($valistep{11}) {
	$sc = "annotate_variation_ABCC.pl -filter -dbtype ljb_mt -score_threshold 0 -buildver $buildver -outfile $outfile $alt_input $dbloc";
	if ($keepline){
		$sc=~s/\-filter/\-regionanno/;
	}
	print STDERR "\nNOTICE: Running step 11 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step12
if ($valistep{12}) {
	$sc = "annotate_variation_ABCC.pl -filter -dbtype ljb_lrt -score_threshold 0 -buildver $buildver -outfile $outfile $alt_input $dbloc";
	if ($keepline){
		$sc=~s/\-filter/\-regionanno/;
	}
	print STDERR "\nNOTICE: Running step 12 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

open (FUNCTION, "$outfile.variant_function") or die "Error: cannot read from variant function file: $!\n";

if ($valistep{1}) {
	open (STEP1, "$outfile.exonic_variant_function") or die "Error: cannot read from exonic variant function file: $!\n";
}

if ($valistep{2}) {
	if ($buildver eq 'hg18') {
		open (STEP2, "$outfile.hg18_phastConsElements44way") or die "Error: cannot read from mce file: $!\n";
	} 
	elsif ($buildver eq 'hg19') {
		open (STEP2, "$outfile.hg19_phastConsElements46way") or die "Error: cannot read from mce file: $!\n";
	}
}

if ($valistep{3}) {
	open (STEP3, "$outfile.${buildver}_genomicSuperDups") or die "Error: cannot read from segdup file: $!\n";
}

if ($buildver eq 'hg19') {
	if ($valistep{4}) {
		open (STEP4, "$outfile.hg19_ALL.sites.2010_11_dropped") or die "Error: cannot read from drop file $outfile.hg19_ALL.sites.2010_11_dropped: $!\n";
	}
} 
else {
	if ($valistep{4}) {
		open (STEP4, "$outfile.hg18_CEU.sites.2010_07_dropped") or die "Error: cannot read from drop file $outfile.hg18_CEU.sites.2010_07_dropped: $!\n";
	}
	if ($valistep{5}) {
		open (STEP5, "$outfile.hg18_YRI.sites.2010_07_dropped") or die "Error: cannot read from drop file $outfile.hg18_YRI.sites.2010_07_dropped: $!\n";
	}
	if ($valistep{6}) {
		open (STEP6, "$outfile.hg18_JPTCHB.sites.2010_07_dropped") or die "Error: cannot read from drop file $outfile.hg18_JPTCHB.sites.2010_07_dropped: $!\n";
	}
}

if ($valistep{7}) {
	open (STEP7, "$outfile.${buildver}_snp${verdbsnp}_dropped") or die "Error: cannot read from snp$verdbsnp drop file: $!\n";
}

if ($valistep{8}) {
	open (STEP8, "$outfile.${buildver}_avsift_dropped") or die "Error: cannot read from avsift drop file: $!\n";
}

if ($valistep{9}) {
	open (STEP9, "$outfile.${buildver}_ljb_pp2_dropped") or die "Error: cannot read from ljb_pp2 drop file: $!\n";
}

if ($valistep{10}) {
	open (STEP10, "$outfile.${buildver}_ljb_phylop_dropped") or die "Error: cannot read from ljb_phylop drop file: $!\n";
}

if ($valistep{11}) {
	open (STEP11, "$outfile.${buildver}_ljb_mt_dropped") or die "Error: cannot read from ljb_mt drop file: $!\n";
}

if ($valistep{12}) {
	open (STEP12, "$outfile.${buildver}_ljb_lrt_dropped") or die "Error: cannot read from ljb_lrt drop file: $!\n";
}

my (@allstep);
for my $i (1 .. 3) {
	$allstep[$i] = [];
}

if ($buildver eq 'hg19') {
	if ($valistep{4}) {
		while (<STEP4>) {
			m/^1000g\w+_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 4 : <$_>\n";
			$allstep[4]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
			$allstep[5]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
			$allstep[6]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
} 
else {
	if ($valistep{4}) {
		while (<STEP4>) {
			m/^1000g\w*_ceu\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 4 : <$_>\n";
			$allstep[4]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}

	if ($valistep{5}) {
		while (<STEP5>) {
			m/^1000g\w*_yri\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 5 : <$_>\n";
			$allstep[5]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	
	if ($valistep{6}) {
		while (<STEP6>) {
			m/^1000g\w*_jptchb\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 6 : <$_>\n";
			$allstep[6]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
}

if ($valistep{7}) {
	while (<STEP7>) {
		m/^snp\d+\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 7 : <$_>\n";
		$allstep[7]->{$2} = $1;
	}
}

if ($valistep{8}) {
	while (<STEP8>) {
		m/^avsift\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 8 : <$_>\n";
		$allstep[8]->{$2} = $1;
	}
}

if ($valistep{9}) {
	while (<STEP9>) {
		m/^ljb_pp2\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 9 : <$_>\n";
		$allstep[9]->{$2} = $1;
	}
}

if ($valistep{10}) {
	while (<STEP10>) {
		m/^ljb_phylop\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 10 : <$_>\n";
		$allstep[10]->{$2} = $1;
	}
}

if ($valistep{11}) {
	while (<STEP11>) {
		m/^ljb_mt\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 11 : <$_>\n";
		$allstep[11]->{$2} = $1;
	}
}

if ($valistep{12}) {
	while (<STEP12>) {
		m/^ljb_lrt\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 12 : <$_>\n";
		$allstep[12]->{$2} = $1;
	}
}

print STDERR "NOTICE: Finished loading filterstep database file\n";

open (OUT, ">$outfile.genome_summary.csv") or die "Error: cannot write to output file: $!\n";
open (OUTE, ">$outfile.exome_summary.csv") or die "Error: cannot write to output file: $!\n";

if ($buildver eq 'hg19') {
	print OUT join (',', qw/Func Gene ExonicFunc AAChange Conserved SegDup 1000G_ALL 1000G_ALL 1000G_ALL/, "dbSNP$verdbsnp", qw/SIFT PolyPhen2 LJB_PhyloP LJB_MutationTaster LJB_LRT Chr Start End Ref Obs Otherinfo/), "\n";
	print OUTE join (',', qw/Func Gene ExonicFunc AAChange Conserved SegDup 1000G_ALL 1000G_ALL 1000G_ALL/, "dbSNP$verdbsnp", qw/SIFT PolyPhen2 LJB_PhyloP LJB_MutationTaster LJB_LRT Chr Start End Ref Obs Otherinfo/), "\n";
} else {
	print OUT join (',', qw/Func Gene ExonicFunc AAChange Conserved SegDup 1000G_CEU 1000G_YRI 1000G_JPTCHB/, "dbSNP$verdbsnp", qw/SIFT PolyPhen2 LJB_PhyloP LJB_MutationTaster LJB_LRT Chr Start End Ref Obs Otherinfo/), "\n";
	print OUTE join (',', qw/Func Gene ExonicFunc AAChange Conserved SegDup 1000G_CEU 1000G_YRI 1000G_JPTCHB/, "dbSNP$verdbsnp", qw/SIFT PolyPhen2 LJB_PhyloP LJB_MutationTaster LJB_LRT Chr Start End Ref Obs Otherinfo/), "\n";
}

while (<FUNCTION>) {
	s/[\r\n]+$//;
	m/^(\S+)\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)(.*)/ or die "Error: invalid record found in annovar outputfile: <$_>\n";
	my ($function, $gene, $varstring, $otherinfo) = ($1, $2, $3, $4||'');
	my $exonic;
	if ($function =~ m/\bsplicing\b/ or $function =~ m/\bexonic\b/) {
		$exonic = 1;
	}
	print OUT qq/"$function","$gene"/;
	$exonic and print OUTE qq/"$function","$gene"/;
	
	if (not @{$allstep[1]}) {
		if ($valistep{1}) {
			if (defined ($_ = <STEP1>)) {
				m/^line\d+\t([^\t]+)\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 2: <$_>\n";
				my ($efun, $aachange, $varstring) = ($1, $2, $3);
				my @aachange = split (/:|,/, $aachange);
				if (@aachange >= 5) {
					push @{$allstep[1]}, $varstring, $efun, "$aachange[1]:$aachange[3]:$aachange[4]";
				} else {
					push @{$allstep[1]}, $varstring, $efun, $aachange;		#aachange could be "UNKNOWN"
				}
			}
		}
	}
	
	if (not @{$allstep[2]}) {
		if ($valistep{2}) {
			if (defined ($_ = <STEP2>)) {
				m/^mce\d+way\tScore=(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 2: <$_>\n";
				push @{$allstep[2]}, $2, $1;
			}
		}
	}
	
	if (not @{$allstep[3]}) {
		if ($valistep{3}) {
			if (defined ($_ = <STEP3>)) {
				m/^segdup\tScore=(\S+);\S+\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 3 : <$_>\n";
				push @{$allstep[3]}, $2, ($1>0.01)?sprintf("%.2f", $1):$1;
			}
		}
	}
	
	for my $i (1 .. 3) {			#starting from step 1 to step 7
		my $curstep = $allstep[$i];
		if (@$curstep and $curstep->[0] eq $varstring) {
			if ($i == 1) {
				print OUT qq/,"$curstep->[1]"/, qq/,"$curstep->[2]"/;
				$exonic and print OUTE qq/,"$curstep->[1]"/, qq/,"$curstep->[2]"/;
			} else {
				print OUT qq/,"$curstep->[1]"/;
				$exonic and print OUTE qq/,"$curstep->[1]"/;
			}
			@$curstep = ();
		}
		else {
			if ($i == 1) {
				print OUT ",,";
				$exonic and print OUTE ",,";
			} else {
				print OUT ",";
				$exonic and print OUTE ",";
			}
		}
	}
	
	for my $i (4 .. 12) {		#step 8 to step 11
		if (defined $allstep[$i]->{$varstring}) {
			print OUT qq/,$allstep[$i]->{$varstring}/;
			$exonic and print OUTE qq/,$allstep[$i]->{$varstring}/;
		} else {
			print OUT ",";
			$exonic and print OUTE ",";
		}
	}
	
	my @varstring = split (/\s+/, $varstring);
	$otherinfo =~ s/^\s+//;
	my @otherinfo = split (/\t/, $otherinfo);
	for my $i (0 .. @otherinfo-1) {
		$otherinfo[$i] = qq/"$otherinfo[$i]"/;
	}

	print OUT ',', join (',', @varstring), ",", join (',', @otherinfo), "\n";
	$exonic and print OUTE ',', join (',', @varstring), ",", join (',', @otherinfo), "\n";
}

print STDERR "NOTICE: Final whole-genome summary was written to $outfile.genome_summary.csv file\n";
print STDERR "NOTICE: Final whole-exome summary was written to $outfile.exome_summary.csv file\n";

if ($remove) {
	unlink ("$outfile.variant_function", "$outfile.exonic_variant_function", "$outfile.hg18_phastConsElements44way", "$outfile.hg19_phastConsElements46way", 
		"$outfile.${buildver}_genomicSuperDups", "$outfile.hg19_ALL.sites.2010_11_dropped", "$outfile.hg18_CEU.sites.2010_07_dropped", "$outfile.hg18_YRI.sites.2010_07_dropped",
		"$outfile.hg18_JPTCHB.sites.2010_07_dropped", "$outfile.${buildver}_snp${verdbsnp}_dropped", "$outfile.${buildver}_avsift_dropped");
}

sub checkFileExistence {
	my @file = ("${buildver}_refGene.txt", "${buildver}_refLink.txt", "${buildver}_refGeneMrna.fa", "${buildver}_genomicSuperDups.txt", 
		"${buildver}_snp$verdbsnp.txt", "${buildver}_avsift.txt", "${buildver}_ljb_pp2.txt", "${buildver}_ljb_phylop.txt", "${buildver}_ljb_mt.txt", "${buildver}_ljb_lrt.txt");
	if ($buildver eq 'hg18') {
		push @file, "${buildver}_phastConsElements44way.txt";
		push @file, "${buildver}_CEU.sites.2010_07.txt", "${buildver}_YRI.sites.2010_07.txt", "${buildver}_JPTCHB.sites.2010_07.txt";
	} elsif ($buildver eq 'hg19') {
		push @file, "${buildver}_phastConsElements46way.txt";
		push @file, "${buildver}_ALL.sites.2010_11.txt";
	}
	for my $i (0 .. @file-1) {
		my $dbfile = File::Spec->catfile ($dbloc, $file[$i]);
		-f $dbfile or die "Error: the required database file $dbfile does not exist. Please download it via -downdb argument by annotate_variation_ABCC.pl.\n";
	}
}

=head1 SYNOPSIS

 summarize_annovar.pl [arguments] <query-file> <database-location>

 Optional arguments:
		-h, --help				print help message
		-m, --man				print complete documentation
		-v, --verbose				use verbose output
		    --outfile <string>		output file name prefix
		    --buildver <string>		genome build version (default: hg18)
		    --remove			remove all temporary files
		    --verdbsnp <int>		dbSNP version to use (default: 130)

 Function: automatically run a pipeline on a list of variants and summarize 
 their functional effects in a comma-delimited file, to be opened by Excel for 
 manual filtering
 
 Example: summarize_annovar.pl ex2.human humandb/
 
 Version: $LastChangedDate: 2011-10-02 21:27:43 -0700 (Sun, 02 Oct 2011) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--outfile>

the prefix of output file names

=item B<--buildver>

specify the genome build version

=item B<--remove>

remove all temporary files. By default, all temporary files will be kept for 
user inspection, but this will easily clutter the directory.

=item B<--checkfile>

check to make sure that database files exist, before executing the current 
program.

=back

=head1 DESCRIPTION

ANNOVAR is a software tool that can be used to functionally annotate a list of 
genetic variants, possibly generated from next-generation sequencing 
experiments. For example, given a whole-genome resequencing data set for a human 
with specific diseases, typically around 3 million SNPs and around half million 
insertions/deletions will be identified. Given this massive amounts of data (and 
candidate disease- causing variants), it is necessary to have a fast algorithm 
that scans the data and identify a prioritized subset of variants that are most 
likely functional for follow-up Sanger sequencing studies and functional assays.


ANNOVAR is freely available to the community for non-commercial use. For 
questions or comments, please contact kai@openbioinformatics.org.


=cut
