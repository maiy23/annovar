#!/usr/bin/perl
#$ver = 10 //arbritarily started at 10
#version10 includes kaviar and hrcr1
#==> hg19_kaviar_20150923.txt <==
#Chr    Start   End     Ref     Alt     Kaviar_AF       Kaviar_AC       Kaviar_AN
# 1       1       1       0       0       .       .       26028
# 1       10001   10001   T       C       3.84e-05        1       26028
# 1       10002   10002   A       AT      3.84e-05        1       26028
# 1       10002   10002   A       C       0.0001153       3       26028
# 1       10002   10002   A       T       3.84e-05        1       26028
# 1       10002   10002   -       T       3.84e-05        1       26028
# 1       10003   10003   A       C       3.84e-05        1       26028
# 1       10003   10003   A       T       7.68e-05        2       26028
# 1       10004   10004   C       A       3.84e-05        1       26028

# ==> hg19_hrcr1.txt <==
# #Chr    Start   End     Ref     Alt     HRC_AF  HRC_AC  HRC_AN  HRC_non1000G_AF HRC_non1000G_AC HRC_non1000G_AN
# 1       13380   13380   C       G       7.69515e-05     5       64976   0       0       59986
# 1       16071   16071   G       A       0.000123122     8       64976   0       0       59986
# 1       16141   16141   C       T       0.000138513     9       64976   6.66822e-05     4       59986
# 1       16280   16280   T       C       0.000661783     43      64976   3.33411e-05     2       59986
# 1       49298   49298   T       C       0.640021        41586   64976   0.636315        38170   59986
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
use Cwd;
#ABCC Global Modifications Here
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use ABCC_moo;
my $bin=`dirname $0`;chomp $bin;
umask(0000);
our (%gene_ori, $relpos, $header, $keepline,%filtermap,%filtererr, $novar,$runPTM , $hasHeader,$protDB);#filtermap keeps track of the order for keepline
#### END ABCC GLobal Modifications
our $VERSION = 			'$Revision: 491 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-10-05 09:10:39 -0700 (Wed, 05 Oct 2011) $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $separate, $batchsize, $dbtype, $neargene, $genomebinsize, $geneanno, $regionanno, $filter, $downdb, $buildver, $score_threshold, $normscore_threshold, $minqueryfrac, $expandbin, $splicing_threshold,
	$maf_threshold, $chromosome, $zerostart, $rawscore, $memfree, $memtotal, $sift_threshold, $gff3dbfile, $genericdbfile, $vcfdbfile, $time, $wget, $precedence,$correctAllele, $force,
	$webfrom, $colsWanted, $comment, $scorecolumn, $transfun, $exonsort, $avcolumn, $bedfile, $hgvs, $reverse, $indexfilter_threshold ,$silent,$extype ,$collapseOnly, $noheader,$match);
our (%valichr, $dbtype1, $headersWanted);
our (@precedence, @colsWanted, @avcolumn);
sub printerr;			#declare a subroutine
our $cwd=`pwd`;chomp $cwd;
our %codon1 = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"*", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codon3 = (TTT=>"Phe", TTC=>"Phe", TCT=>"Ser", TCC=>"Ser", TAT=>"Tyr", TAC=>"Tyr", TGT=>"Cys", TGC=>"Cys", TTA=>"Leu", TCA=>"Ser", TAA=>"*", TGA=>"*", TTG=>"Leu", TCG=>"Ser", TAG=>"*", TGG=>"Trp", CTT=>"Leu", CTC=>"Leu", CCT=>"Pro", CCC=>"Pro", CAT=>"His", CAC=>"His", CGT=>"Arg", CGC=>"Arg", CTA=>"Leu", CTG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", ATT=>"Ile", ATC=>"Ile", ACT=>"Thr", ACC=>"Thr", AAT=>"Asn", AAC=>"Asn", AGT=>"Ser", AGC=>"Ser", ATA=>"Ile", ACA=>"Thr", AAA=>"Lys", AGA=>"Arg", ATG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"Arg", GTT=>"Val", GTC=>"Val", GCT=>"Ala", GCC=>"Ala", GAT=>"Asp", GAC=>"Asp", GGT=>"Gly", GGC=>"Gly", GTA=>"Val", GTG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonfull = (TTT=>"Phenylalanine", TTC=>"Phenylalanine", TCT=>"Serine", TCC=>"Serine", TAT=>"Tyrosine", TAC=>"Tyrosine", TGT=>"Cysteine", TGC=>"Cysteine", TTA=>"Leucine", TCA=>"Serine", TAA=>"Stop", TGA=>"Stop", TTG=>"Leucine", TCG=>"Serine", TAG=>"Stop", TGG=>"Tryptophan", CTT=>"Leucine", CTC=>"Leucine", CCT=>"Proline", CCC=>"Proline", CAT=>"Histidine", CAC=>"Histidine", CGT=>"Arginine", CGC=>"Arginine", CTA=>"Leucine", CTG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", ATT=>"Isoleucine", ATC=>"Isoleucine", ACT=>"Threonine", ACC=>"Threonine", AAT=>"Asparagine", AAC=>"Asparagine", AGT=>"Serine", AGC=>"Serine", ATA=>"Isoleucine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Arginine", ATG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Arginine", GTT=>"Valine", GTC=>"Valine", GCT=>"Alanine", GCC=>"Alanine", GAT=>"Aspartic acid", GAC=>"Aspartic acid", GGT=>"Glycine", GGC=>"Glycine", GTA=>"Valine", GTG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine");
our %codonr1 = (UUU=>"F", UUC=>"F", UCU=>"S", UCC=>"S", UAU=>"Y", UAC=>"Y", UGU=>"C", UGC=>"C", UUA=>"L", UCA=>"S", UAA=>"*", UGA=>"*", UUG=>"L", UCG=>"S", UAG=>"*", UGG=>"W", CUU=>"L", CUC=>"L", CCU=>"P", CCC=>"P", CAU=>"H", CAC=>"H", CGU=>"R", CGC=>"R", CUA=>"L", CUG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", AUU=>"I", AUC=>"I", ACU=>"T", ACC=>"T", AAU=>"N", AAC=>"N", AGU=>"S", AGC=>"S", AUA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", AUG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GUU=>"V", GUC=>"V", GCU=>"A", GCC=>"A", GAU=>"D", GAC=>"D", GGU=>"G", GGC=>"G", GUA=>"V", GUG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codonr3 = (UUU=>"Phe", UUC=>"Phe", UCU=>"Ser", UCC=>"Ser", UAU=>"Tyr", UAC=>"Tyr", UGU=>"Cys", UGC=>"Cys", UUA=>"Leu", UCA=>"Ser", UAA=>"*", UGA=>"*", UUG=>"Leu", UCG=>"Ser", UAG=>"*", UGG=>"Trp", CUU=>"Leu", CUC=>"Leu", CCU=>"Pro", CCC=>"Pro", CAU=>"His", CAC=>"His", CGU=>"Arg", CGC=>"Arg", CUA=>"Leu", CUG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", AUU=>"Ile", AUC=>"Ile", ACU=>"Thr", ACC=>"Thr", AAU=>"Asn", AAC=>"Asn", AGU=>"Ser", AGC=>"Ser", AUA=>"Ile", ACA=>"Thr", AAA=>"Lys", AGA=>"Arg", AUG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"Arg", GUU=>"Val", GUC=>"Val", GCU=>"Ala", GCC=>"Ala", GAU=>"Asp", GAC=>"Asp", GGU=>"Gly", GGC=>"Gly", GUA=>"Val", GUG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonrfull = (UUU=>"Phenylalanine", UUC=>"Phenylalanine", UCU=>"Serine", UCC=>"Serine", UAU=>"Tyrosine", UAC=>"Tyrosine", UGU=>"Cysteine", UGC=>"Cysteine", UUA=>"Leucine", UCA=>"Serine", UAA=>"Stop", UGA=>"Stop", UUG=>"Leucine", UCG=>"Serine", UAG=>"Stop", UGG=>"Tryptophan", CUU=>"Leucine", CUC=>"Leucine", CCU=>"Proline", CCC=>"Proline", CAU=>"Histidine", CAC=>"Histidine", CGU=>"Arginine", CGC=>"Arginine", CUA=>"Leucine", CUG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", AUU=>"Isoleucine", AUC=>"Isoleucine", ACU=>"Threonine", ACC=>"Threonine", AAU=>"Asparagine", AAC=>"Asparagine", AGU=>"Serine", AGC=>"Serine", AUA=>"Isoleucine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Arginine", AUG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Arginine", GUU=>"Valine", GUC=>"Valine", GCU=>"Alanine", GCC=>"Alanine", GAU=>"Aspartic acid", GAC=>"Aspartic acid", GGU=>"Glycine", GGC=>"Glycine", GUA=>"Valine", GUG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine");
our %iupac = (R=>'AG', Y=>'CT', S=>'GC', W=>'AT', K=>'GT', M=>'AC', A=>'AA', C=>'CC', G=>'GG', T=>'TT', B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT', '.'=>'-', '-'=>'-');

processArguments ();			#process program arguments, set up default values, check for errors, check for existence of db files

##note to self: if all variants 'fail' (no hit), then check that the inputs are stripped of 'chr' in chromosome positions in each of these subroutines
if ($geneanno) {
	if (!$collapseOnly){annotateQueryByGene ();	}#generate gene-based annoations (classify variants into intergenic, introgenic, non-synonymous, synonymous, UTR, frameshift, etc)
	if ($keepline){
		_collapseANNOVAR("$outfile.variant_function",0);
		#run Protein databases specified in the protDB option
		$protDB.=',PTM' if ($protDB!~/PTM/ && $runPTM);
		if ($protDB && $dbtype1!~/ensGene/){
			my $addon;
			my @protDBs=split(",",$protDB);
			if ($noheader){
				$addon=" -h" ;
			}
			printerr ("perl $bin/protein_impactsIP.pl -I $outfile.variant_function.collapsed -k $addon -D $protDB -v $buildver\n");
			system ("perl $bin/protein_impactsIP.pl -I $outfile.variant_function.collapsed -k $addon -D $protDB -v $buildver\n");
			# foreach my $protdb (@protDBs){
			# 	next if (!$protdb);
			# 	if ($protdb=~/PTM/ ){
			# 		printerr "perl $bin/protein_impacts.pl -I $outfile.variant_function.collapsed -k -o $outfile.$buildver\_PTM $addon\n";next;
			# 		exec ("perl $bin/protein_impacts.pl -I $outfile.variant_function.collapsed -k -o $outfile.$buildver\_PTM $addon\n");
			# 	}else{
			# 		printerr "Running...perl $bin/protein_impacts.pl -I $outfile.variant_function.collapsed -k -o $outfile.$buildver\_$protdb $addon -D $buildver\_$protdb.txt\n";next;
			# 		exec ("perl $bin/protein_impacts.pl -I $outfile.variant_function.collapsed -k -o $outfile.$buildver\_$protdb $addon -D $buildver\_$protdb.txt\n");
			# 	}
			# }
		}else{
			printerr "No protDBs requested\n$protDB\n";
		}
		printerr "NOTICE: Collapsed $outfile.variant_function and $outfile.exonic_variant_function to $outfile.variant_function.collapsed\n";
	}
	printerr "\nANNOVAR Completed:\n\t". `date`;
} elsif ($regionanno) {
	if ($dbtype=~/(encode|rptMask)/i){
		# print "shortcut\n";
		annotateQueryByRegionWithIdx();
	}else{
		annotateQueryByRegion();	#generate region-based annotations (most conserved elements, transcription factor binding sites, etc)
	}
} elsif ($filter ) {
	my $id=filterQuery ();			#generate filter-based annotations (identify variants not reported in variation databases)
	if ($keepline){
		_sortFilterFile($id);
		printerr "NOTICE: Collapsed $id\_filtered and $id\_dropped to $id\n";
	}
	printerr "\nANNOVAR Completed:\n\t". `date`;
} elsif ($downdb) {
	downloadDB ();			#download annotation databases from Internet
}

sub processArguments {
	my @command_line = @ARGV;		#command line argument
	GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'separate'=>\$separate,
	'batchsize=s'=>\$batchsize, 'dbtype=s'=>\$dbtype, 'neargene=i'=>\$neargene, 'genomebinsize=s'=>\$genomebinsize,
	'geneanno'=>\$geneanno, 'regionanno'=>\$regionanno, , 'filter'=>\$filter, 'downdb'=>\$downdb, 'buildver=s'=>\$buildver, 'score_threshold=f'=>\$score_threshold, 
	'normscore_threshold=i'=>\$normscore_threshold,	'minqueryfrac=f'=>\$minqueryfrac, 'expandbin=i'=>\$expandbin, 'splicing_threshold=i'=>\$splicing_threshold,
	'maf_threshold=f'=>\$maf_threshold, 'chromosome=s'=>\$chromosome, 'zerostart'=>\$zerostart, 'rawscore'=>\$rawscore, 'memfree=i'=>\$memfree, 
	'memtotal=i'=>\$memtotal, 'sift_threshold=f'=>\$sift_threshold, 'gff3dbfile=s'=>\$gff3dbfile, 'genericdbfile=s'=>\$genericdbfile, 'vcfdbfile=s'=>\$vcfdbfile,
	'time'=>\$time, 'wget!'=>\$wget, 'precedence=s'=>\$precedence, 'webfrom=s'=>\$webfrom, 'colsWanted=s'=>\$colsWanted, 'comment'=>\$comment,
	'scorecolumn=i'=>\$scorecolumn, 'transcript_function'=>\$transfun, 'exonsort'=>\$exonsort, 'avcolumn=s'=>\$avcolumn, 'bedfile=s'=>\$bedfile, 
	'relpos'=>\$relpos, 'header'=>\$header, 'keepline'=>\$keepline, 'silent'=>\$silent, 'novar'=>\$novar, 'extype'=>\$extype,'collapseOnly'=>\$collapseOnly, 'runPTM'=>\$runPTM, 'noheader'=>\$noheader, 'match'=>\$match, 'correctAllele'=>\$correctAllele, 'force'=>\$force, 'protDB=s'=>\$protDB, 'headersWanted=s'=>\$headersWanted, 
	'hgvs'=>\$hgvs, 'reverse'=>\$reverse, 'indexfilter_threshold=f'=>\$indexfilter_threshold) or pod2usage ();#hv added relpos and header and keepline, silent and novar,$headersWanted 
	
	$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
	$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
	@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
	@ARGV == 2 or pod2usage ("Syntax error:$#ARGV". join (" ",@ARGV));
	$novar ||= '';
	($queryfile, $dbloc) = @ARGV;

	$dbloc =~ s/[\\\/]$//;			#delete the trailing / or \ sign as part of the directory name
	if (defined $batchsize) {
		$batchsize =~ s/k$/000/;
		$batchsize =~ s/m$/000000/;
		$batchsize =~ m/^\d+$/ or pod2usage ("Error: the --batchsize argument must be a positive integer (suffix of k or m is okay)");
	} else {
		$batchsize = 1_000_000;
	}
	if (defined $genomebinsize) {
		$genomebinsize =~ s/k$/000/;
		$genomebinsize =~ s/m$/000000/;
		$genomebinsize =~ m/^\d+$/ or pod2usage ("Error: the --genomebinsize argument must be a positive integer (suffix of k or m is okay)");
		$genomebinsize > 1000 or pod2suage ("Error: the --genomebinsize argument must be larger than 1000");
	} else {
		if ($geneanno) {
			$genomebinsize = 100_000;		#gene usually span large genomic regions
		} else {
			$genomebinsize = 10_000;		#MCE, TFBS, miRNA, etc are small genomic regions
		}
	}

	$verbose ||= 0;			#when it is not specified, it is zero
	$neargene ||= 1_000;		#for upstream/downstream annotation of variants, specify the distance threshold between variants and genes
	$expandbin ||= int(2_000_000/$genomebinsize);		#for gene-based annotations, when intergenic variants are found, expand to specified number of nearby bins to find closest genes
	$outfile ||= $queryfile;	#specify the prefix of output file names

	#set up log file
	if ($downdb) {
		if (not -d $dbloc) {
			mkdir ($dbloc) or die "Error: the directory $dbloc does not exist and cannot be created\n";
		}
		my $errfile = File::Spec->catfile ($dbloc, "annovar_downdb.log");
		open (LOG, ">$errfile") or die "Error: cannot write LOG information to log file $errfile: $!\n";
	} else {
		my $dbaddon='';##hv added next block 8/17/15 so log files will not be overwritten
		if ($genericdbfile){
			$dbaddon="_$genericdbfile";
			$dbaddon=~s/(\.txt|hg19_|\.gff|\.text|\.csv)//g;
		}elsif ($dbtype){
			$dbaddon="_$dbtype";
		}elsif($geneanno){
			$dbaddon="_geneanno";
		}
		open (LOG, "| tee $outfile$dbaddon.log") or die __LINE__. "REINSTATE!!! -> Error: cannot write LOG information to log file $outfile.log: $!\n";
		# print LOG ("\t\t\t\t\t--------->Reinstate this line " . __LINE__."\n");
	}
	printerr "ANNOVAR Running on host($novar): ". `hostname`;
	printerr "ANNOVAR Version:\n\t", q/$LastChangedDate: 2011-10-05 09:10:39 -0700 (Wed, 05 Oct 2011) $/, "\n";
	printerr "ANNOVAR Information:\n\tFor questions, comments, documentation, bug reports and program update, please visit http://www.openbioinformatics.org/annovar/\n";
	printerr "ANNOVAR Command:\n\t$0 @command_line\n";
	printerr "ANNOVAR Started:\n\t". `date`. "\t(with modified hash annots)\n";
	
	my $num = 0;
	$geneanno and $num++;
	$downdb and $num++;
	$filter and $num++;
	$regionanno and $num++;
	$num <= 1 or pod2usage ("Error in argument: please specify only one of --geneanno, -regionanno, --downdb, --filter");
	if (not $num) {
		$geneanno++;
		printerr "NOTICE: The --geneanno operation is set to ON by default\n";
	}
	
	my %dbtype1 = ('gene'=>'refGene', 'refgene'=>'refGene', 'knowngene'=>'knownGene', 'ensgene'=>'ensGene', 'band'=>'cytoBand', 'cytoband'=>'cytoBand', 'tfbs'=>'tfbsConsSites', 'mirna'=>'wgRna',
			'mirnatarget'=>'targetScanS', 'segdup'=>'genomicSuperDups', 'omimgene'=>'omimGene', 'gwasCatalog'=>'gwasCatalog', 
			'1000g_ceu'=>'CEU.sites.2009_04', '1000g_yri'=>'YRI.sites.2009_04', '1000g_jptchb'=>'JPTCHB.sites.2009_04', 
			'1000g2010_ceu'=>'CEU.sites.2010_03', '1000g2010_yri'=>'YRI.sites.2010_03', '1000g2010_jptchb'=>'JPTCHB.sites.2010_03',
			'1000g2010jul_ceu'=>'CEU.sites.2010_07', '1000g2010jul_yri'=>'YRI.sites.2010_07', '1000g2010jul_jptchb'=>'JPTCHB.sites.2010_07',
			'1000g2010nov_all'=>'ALL.sites.2010_11', '1000g2011may_all'=>'ALL.sites.2011_05','ALL_sites_2011_05'=>'ALL.sites.2011_05'
			);
	
	if ($geneanno) {
		$dbtype ||= 'refGene';
		$dbtype1 = $dbtype1{$dbtype} || $dbtype;
		#$dbtype1 =~ m/^(refGene|knownGene|ensGene)$/ or pod2usage ("Error: the gene-based annotation procedure currently only support -dbtype of refGene, knownGene and ensGene");	#commented 2011feb18
	} elsif ($regionanno) {
		defined $dbtype or pod2usage ("Error in argument: please specify --dbtype (required for the --regionanno operation)");
		$dbtype1 = $dbtype1{$dbtype} || $dbtype;
		if ($dbtype =~ m/^mce(\d+)way/) {			#added 2010Feb16
			$dbtype1 = "phastConsElements$1way";
		}
		if ($dbtype1 eq 'gff3') {
			defined $gff3dbfile or pod2usage ("Error in argument: please specify --gff3dbfile for the --dbtype of 'gff3'");
		}
	} elsif ($filter) {
		defined $dbtype  or die ("Error in argument: please specify --dbtype (required for the --filter operation) ($dbtype)");
		# if ($dbtype =~ m/(ALL|AFR|AMR|EAS|EUR|SAS)\.sites|exac0\d|^1000gALL.sites.[\.\d\-\_]*|pp2|cosmicdb|hapmap|avsift|generic|1000g_(ceu|yri|jptchb)|rbpvar|1000g2010_(ceu|yri|jptchb)|nci60|1000g20\d\d[a-z]{3}_[a-z]+|regulome\d+|snp\d+|vcf|icgc|(ljb\d{0,2}_\w+|hom_only|cg69|clinvar|cosmic70|ExAC|viz\d+)|cadd|esp6500|ESP6500|siftv63|Provean|hrcr|kaviar_|dbsnp_v\d+$/ and !$force){
		# 		pod2usage ("Error in argument: the specified --dbtype $dbtype is not valid for --filter operation (valid ones are '1000g_ceu', '1000g2010_yri', 'snp129', 'avsift', 'vcf', 'generic', etc)");
		# }
		$dbtype1 = $dbtype1{$dbtype} || $dbtype;
		if ($dbtype1 eq 'generic') {

			defined $genericdbfile or pod2usage ("Error in argument: please specify --genericdbfile for the --dbtype of 'generic'");
		}
		if ($dbtype eq 'vcf') {
			defined $vcfdbfile or pod2usage ("Error in argument: please specify --vcfdbfile for the --dbtype of 'vcf'");
		}
	} elsif ($downdb) {
		defined $dbtype and pod2usage ("Error in argument: please do not specify --dbtype for the --downdb operation");
		$dbtype1 = $dbtype1{$queryfile} || $queryfile;
	}
	if ($match && !$regionanno){
		pod2usage ("Error in argument: --match is only supported with the --regionanno operation\n");
	}
	if (not $buildver) {
		$buildver = 'hg19';
		printerr "NOTICE: The --buildver is set as 'hg19' by default\n";
	}
	
	if (defined $score_threshold) {
		$score_threshold >= 0 or pod2usage ("Error in argument: the --score_threshold must be a positive number or zero (you specified $score_threshold)");
		$geneanno || $downdb and pod2usage ("Error in argument: the --score_threshold is not useful for --geneanno or --downdb operations");
	}
	if ($normscore_threshold) {
		$normscore_threshold <= 1000 or pod2usage ("Error in argument: the --normscore_threshold must be between 0 and 1000 (you specified $normscore_threshold)");
		$regionanno or pod2usage ("Error in argument: the --score_threshold is supported only for the --regionanno operation");
	}
	
	if ($zerostart) {
		printerr ("WARNING: -zerostart is OBSELETE!!! Do NOT use it!!! Try to reformat your input!!!\n");
		printerr ("WARNING: -zerostart will be permenantly removed in future ANNOVAR release, as many users mistakenly treat it as UCSC's half-open zero-start format.\n");
	}
	
	if (defined $sift_threshold) {
		$filter or pod2usage ("Error in argument: the --sift_threshold is supported only for the --filter operation");
		$dbtype1 eq 'avsift' or pod2usage ("Error in argument: the --sift_threshold argument can be used only if '--dbtype avsift' is used");
		$sift_threshold >= 0 and $sift_threshold <= 1 or pod2usage ("Error in argument: the --sift_threshold must be between 0 and 1 inclusive");
	} else {
		$dbtype1 eq 'avsift' and printerr "NOTICE: The --sift_threshold is set as 0.05 by default\n";
		$sift_threshold = 0.05;
	}
	
	if (defined $indexfilter_threshold) {
		$filter or pod2usage ("Error in argument: the --indexfilter_threshold is supported only for the --filter operation");
		$indexfilter_threshold >= 0 and $indexfilter_threshold <= 1 or pod2usage ("Error in argument: the --indexfilter_threshold must be between 0 and 1 inclusive");
	} else {
		$indexfilter_threshold = 0.5;
	}
	
	#operation-specific argument
	if (defined $splicing_threshold) {
		$geneanno or pod2usage ("Error in argument: the --splicing_threshold is supported only for the --geneanno operation");
	} else {
		$splicing_threshold = 2;	#for splicing annotation, specify the distance threshold between variants and exon/intron boundaries
	}
	if (defined $maf_threshold) {
		$filter or pod2usage ("Error in argument: the --maf_threshold is supported only for the --filter operation");
	} else {
		$maf_threshold = 0;		#for filter-based annotations on 1000 Genomes Project data, specify the MAF threshold to be used in filtering
	}
	if (defined $minqueryfrac) {
		$regionanno or pod2usage ("Error in argument: the --minqueryfrac is supported only for the --regionanno operation");
	} else {
		$minqueryfrac = 0;		#minimum query overlap to declare a "match" with database records
	}
	if (defined $gff3dbfile) {
		$dbtype eq 'gff3' or pod2usage ("Error in argument: the --gff3dbfile argument can be used only if '--dbtype gff3' is used");
		$geneanno or $regionanno or pod2usage ("Error in argument: the --gff3dbfile argument is supported only for the --geneanno or --regionanno operation");
	}
	if (defined $bedfile) {
		$dbtype eq 'bed' or pod2usage ("Error in argument: the --bedfile argument can be used only if '--dbtype bed' is used");
		$regionanno or pod2usage ("Error in argument: the --bedfile argument is supported only for the --regionanno operation");
	}
	if (defined $genericdbfile) {
		$filter or pod2usage ("Error in argument: the --genericdbfile argument is supported only for the --filter operation") ;
	}
	if (defined $wget) {
		$downdb or pod2usage ("Error in argument: the --wget argument is supported only for the --downdb operation");
	} else {
		$wget = 1;			#by default, use wget for downloading files from Internet
	}
	if (defined $precedence) {
		$geneanno or pod2usage ("Error in argument: the --precedence argument is supported only for the --geneanno operation");
		@precedence = split (/,/, $precedence);
		@precedence >= 2 or pod2usage ("Error in argument: the --precedence argument should be comma delimited");
		for my $i (0 .. @precedence-1) {
			$precedence[$i] =~ m/^(exonic|intronic|splicing|utr5|utr3|upstream|downstream|splicing|ncrna)$/ or pod2usage ("Error in argument: the --precedence argument contains invalid keywords (valid ones are exonic|intronic|splicing|utr5|utr3|upstream|downstream|splicing)");
		}
	}
	
	if (defined $colsWanted) {
		($regionanno or ($filter && $dbtype=~/exac\d+/) or ($filter && $force)) or pod2usage ("Error in argument: the --colWanted argument is supported only for the --regionanno operation");
		if (lc $colsWanted eq 'all') {
			@colsWanted = ('all');
		} elsif (lc $colsWanted eq 'none') {
			@colsWanted = ('none');
		} else {
			@colsWanted = split (/,/, $colsWanted);
			for my $i (0 .. @colsWanted-1) {
				$colsWanted[$i]=~m/^\d+$/ or pod2usage ("Error in argument: the --colsWanted argument ($colsWanted) must be a list of comma delimited numbers or be 'all' or be 'none'");
			}
		}
	}
	
	if (defined $scorecolumn) {
		$regionanno or pod2usage ("Error in argument: the --scorecolumn argument is supported only for the --regionanno operation");
	}
	
	if ($exonsort) {
		$geneanno or pod2usage ("Error in argument: the --exonsort argument is supported only for the --geneanno operation");
	}
	
	if (defined $avcolumn) {
		$avcolumn =~ m/^\d+,\d+,\d+,\d+,\d+$/ or pod2usage ("Error in argument: the --avcolumn argument must be five integer numbers separated by comma");
		@avcolumn = split (/,/, $avcolumn);
		@avcolumn = map {$_-1} @avcolumn;
	} else {
		@avcolumn = (0..4);		#by default, the first five columns are the required AVINPUT information
	}
	
	if (defined $webfrom) {
		if ($webfrom ne 'ucsc' and $webfrom ne 'annovar') {
			$webfrom =~ m#^(http://|ftp://)# or pod2usage ("Error: the --webfrom argument needs to be 'ucsc', 'annovar', or a URL");
		}
	}
	
	$maf_threshold >= 0 and $maf_threshold <= 0.5 or pod2usage ("Error in argument: the --maf_threshold must be between 0 and 0.5 (you specified $maf_threshold)");
	$minqueryfrac >= 0 and $minqueryfrac <= 1 or pod2usage ("Error in argument: the --minqueryfrac must be between 0 and 1 (you specified $minqueryfrac)");
	$memfree and $memfree >= 100_000 || pod2usage ("Error in argument: the --memfree argument must be at least 100000 (in the order of kilobytes)");
	$memtotal and $memtotal >= 100_000 || pod2usage ("Error in argument: the --memtotal argument must be at least 100000 (in the order of kilobytes)");
	
	if ($chromosome) {
		my @chr = split (/,/, $chromosome);
		for my $i (0 .. @chr-1) {
			if ($chr[$i] =~ m/^(\d+)-(\d+)$/) {
				for my $j ($1 .. $2) {
					$valichr{$j}++;
				}
			} else {
				$valichr{$chr[$i]}++;
			}
		}
		printerr "NOTICE: These chromosomes in database will be examined: ", join (",", sort keys %valichr), "\n";
	}
}


sub annotateQueryByGene {
	my ($queryfh);							#query file handle
	my ($totalquerycount, $totalinvalidcount, $batchcount) = qw/0 0 1/;
	open ($queryfh, $queryfile) or die "Error: cannot read from --queryfile ($queryfile): $!\n";
	
	open (OUT, ">$outfile.variant_function") or die "Error: cannot write to output file $outfile.variant_function: $!\n";
#	if ($keepline){
#		open (EXONIC, '>&OUT') or die "Cannot append to filehandle OUT at line".__LINE__."\n";
#	}else{
		open (EXONIC, ">$outfile.exonic_variant_function") or die "Error: cannot write to output file $outfile.exonic_variant_function: $!\n" if (!$novar);#cannot calculate a protein impact if there is no vars!
#	}
	open (NEW,">$outfile.valid_input") or die "Error cannot write to $outfile.valid_input:$!\n";
	open (INVALID, ">$outfile.invalid_input") or die "Error: cannot write to output file $outfile.invalid_input: $!\n";
	open (INTRON,">$outfile.notexon.input") or die "Error: cannot write to output file $outfile.notexon.input: $!\n" ;

	my ($genedb, $geneidmap, $cdslen, $mrnalen) =readUCSCGeneAnnotation ($dbloc);

	$time and printerr "NOTICE: Current time (before examining variants) is ", scalar (localtime), "\n";
	while (1) {
		my ($linecount, $invalidcount);
		#hv added next block
		my $name_type="Gene";
		if ($transfun){$name_type="transcript";}
		#try to find out if input has has headers 
		$hasHeader=`head -n1 $queryfile | grep -e '^#' `;
		# if (!$hasHeader){
		# 	$hasHeader=`head -n1 $queryfile | grep -vP '^(chr){0,}[\\dXYMT]{1,2}\\t\\d+\\t\\d+'`;
		# 	$hasHeader="#".$hasHeader if ($hasHeader);
		# }
		if (!$relpos){
			($linecount,$invalidcount)= newprocessNextQueryBatchByGene ($queryfh, $batchsize, $genedb, $geneidmap, $cdslen, $mrnalen);
		}else{
			($linecount,$invalidcount) = _newprocessNextQueryBatchByGene_withExtras ($queryfh, $batchsize, $genedb, $geneidmap, $cdslen, $mrnalen, $separate);
		}
		#end block = newprocessNextQueryBatchByGene ($queryfh, $batchsize, $genedb, $geneidmap, $cdslen, $mrnalen);
		$totalquerycount += $linecount;
		$totalinvalidcount += $invalidcount;
		$linecount == $batchsize or last;
		$batchcount++;
		printerr "NOTICE: Begin processing batch $batchcount (each batch contains $batchsize variants)\n" ;
	}
	close (INVALID);
	close (EXONIC);
	close NEW;
	close (OUT);
	close (INTRON);
	close ($queryfh);
	$time and printerr "NOTICE: Current time (after examining variants) is ", scalar (localtime), "\n";

	$totalinvalidcount or unlink ("$outfile.invalid_input");	#delete the file as it is empty
	printerr "NOTICE: Finished gene-based annotation on $totalquerycount genetic variants in $queryfile";
	$totalinvalidcount and printerr " (including $totalinvalidcount with invalid format written to $outfile.invalid_input)";
	!$totalinvalidcount and printerr " (no invalid entries)";
	printerr "\n";
	printerr "NOTICE: Output files were written to $outfile.variant_function, $outfile.exonic_variant_function\n";
	printerr "ANNOVAR Completed:\n\t".scalar (localtime)."\n";
}

sub newprocessNextQueryBatchByGene {
	my ($queryfh, $batchsize, $genedb, $geneidmap, $cdslen, $mrnalen) = @_;
	my (%refseqvar);
	
	my ($chr, $start, $end, $ref, $obs);
	my ($name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend, $name2); 
	my ($invalid);
	my ($linecount, $invalidcount) = qw/0 0/;
	
	for my $i (1 .. $batchsize) {					#process up to batchsize variants
		my $nextline = <$queryfh>;				#read the next line in variant file
		defined $nextline or last;
		$nextline =~ s/[\r\n]+$//;
			$linecount++;	#linecount does not include the comment line
		if ($nextline =~ m/^#/ and $comment) {			#comment line start with #, do not include this is $linecount
			print OUT "#comment\t$nextline\n";
			next;
		}elsif ($nextline=~/^#/){
			next;
		}elsif ($i==1 && $hasHeader){
			print OUT "$hasHeader";
			next;
		}
	
							
		$invalid = 0;
		
		my @nextline = split (/\s+/, $nextline);
		($chr, $start, $end, $ref, $obs) = @nextline[@avcolumn];
		if ( not (defined $chr and defined $start and defined $end and defined $ref and defined $obs)) {
			$invalid++;
		} else {
			($ref, $obs) = (uc $ref, uc $obs);
			$zerostart and $start++;
			$chr =~ s/^chr//i;
			if ($chr =~ m/[^\w]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {
				$invalid++;
			} elsif ($ref eq '-' and $obs eq '-' 		#both are empty allele
				or $ref =~ m/[^ACTG0\-]/ 		#non-standard nucleotide code
				or $obs =~ m/[^ACGT0\-]/ 		#non-standard nucleotide code
				or $start =~ m/[^\d]/ 			#start is not a number
				or $end =~ m/[^\d]/ 			#end is not a number
				or $start > $end			#start is more than end
				or $ref ne '0' and $end-$start+1 != length ($ref) 	#length mismatch with ref
				or $ref eq '-' and $start != $end	#length mismatch for insertion
				) {
				$invalid++;

			}
		}



		if ($invalid) {
			print INVALID $nextline, "\n";			#invalid record found
			print OUT "ERR\t\t\t\t$nextline\n" if ($keepline);
			$invalidcount++;
			next;
		}elsif ($keepline){
			print NEW "$nextline\n";
		}
		
		my (%intronic, %utr5, %utr3, %exonic, %upstream, %downstream, %ncrna, %intergenic, %splicing, %splicing_anno);
		my $foundgenic;						#variant found in genic region (between start and end position of a gene in genome)
		my ($distl, $distr, $genel, $gener);			#for intergenic variant, the distance and gene name to the left and right side of gene
		my $bin1 = int ($start/$genomebinsize)-1;		#start bin
		$bin1 < 0 and $bin1=0;
		my $bin2 = int ($end/$genomebinsize)+1;			#end bin (usually same as start bin, unless the query is really big that spans multiple megabases)
		
		while (not exists $genedb->{$chr, $bin1} and $bin1 > int ($start/$genomebinsize)-$expandbin) {		#examine at least 5 bins (by default 5Mb) to the left to make sure that a gene is found in the bin
			$bin1 > 0 or last;
			$bin1--;
		}
		
		while (not exists $genedb->{$chr, $bin2} and $bin2 < int ($end/$genomebinsize)+$expandbin) {		#examine at least 5 bins (by default 5Mb) to the right to make sure that a gene is found in the bin
			$bin2++;
		}

		my (%seen);
		for my $nextbin ($bin1 .. $bin2) {
			exists $genedb->{$chr, $nextbin} or next;		#this genome bin has no annotated gene (a complete intergenic region)
			for my $nextgene (@{$genedb->{$chr, $nextbin}}) {	#when $genedb->{$chr, $nextbin} is undefined, this automatically create an array!!!
				($name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend, $name2) = @$nextgene;
				defined $name2 or printerr "WARNING: name2 field is not provided for transcript $name (start=$txstart end=$txend)\n" and $name2='';
				$seen{$name, $txstart} and next;		#name and txstart uniquely identify a transcript and chromosome position (sometimes same transcript may map to two nearby positions, such as nearby segmental duplications)
				$seen{$name, $txstart}++;			#a transcript may be in two adjacent bins, so if one is already scanned, there is no need to work on it again
				my $current_ncRNA;
				
				if ($transfun) {	#variant_function output contains transcript name, rather than gene name
					$name2 = $name;
				}
				
				if (not $foundgenic) {				#this variant has not hit a genic region yet
					if ($start > $txend) {
						defined $distl or $distl = $start-$txend and $genel=$name2;
						$distl > $start-$txend and $distl = $start-$txend and $genel=$name2;	#identify left closest gene
					}
		
					if ($end < $txstart) {
						defined $distr or $distr = $txstart-$end and $gener=$name2;
						$distr > $txstart-$end and $distr = $txstart-$end and $gener=$name2;	#identify right closest gene
					}
				}
				
				if ($end < $txstart) {
					#query ---
					#gene		<-*----*->
					$foundgenic and last;					#if found a genic annotation already, end the search of the bins
					if ($end > $txstart - $neargene) {
						if ($dbstrand eq '+') {
							$upstream{$name2}++;
						} else {
							$downstream{$name2}++;
						}
					} else {
						last;						#if transcript is too far away from end, end the search of the bins
					}
				} elsif ($start > $txend) {
					#query            ---
					#gene  <-*----*->
					if (not $foundgenic and $start < $txend + $neargene) {
						if ($dbstrand eq '+') {
							$downstream{$name2}++;
						} else {
							$upstream{$name2}++;
						}
					}
				#} elsif ($cdsstart == $cdsend+1) {				#non-coding RNA (could be microRNA, or could be due to lack of CDS annotation for mRNA such as NR_026730 or BC039000). Previously we already did cdsstart++ so here the cdsstart is more than cdsend
				#	if ($start >= $txstart and $start <= $txend or $end >= $txstart and $end <= $txend or $start <= $txstart and $end >= $txend) {
				#		$ncrna{$name2}++;
				#		$foundgenic++;
				#	}
				} else {							#query overlaps with coding region of gene
					
					#change in 2011jul24: handle ncRNA and protein coding gene together but with the ncRNA flag when printing out results in the future
					if ($cdsstart == $cdsend+1) {				#non-coding RNA (could be microRNA, or could be due to lack of CDS annotation for mRNA such as NR_026730 or BC039000). Previously we already did cdsstart++ so here the cdsstart is more than cdsend
						if ($start >= $txstart and $start <= $txend or $end >= $txstart and $end <= $txend or $start <= $txstart and $end >= $txend) {
							$ncrna{$name2}++;
							$foundgenic++;
						}
						
						#now treat this ncRNA as if it is a protein-coding gene
						($cdsstart, $cdsend) = ($txstart, $txend);
						$current_ncRNA++;		#current transcript is a noncoding transcript
					}
					
					
					
					my ($lenintron, $lenexon) = (0, 0);			#cumulative intron/exon length at a given exon
					my ($rcdsstart, $rvarstart, $rvarend);			#start of coding and variant in reference mRNA sequence
					my @exonstart = @$exonstart;
					my @exonend = @$exonend;
					my $foundexonic;
					if ($dbstrand eq '+') {					#forward strand, search from left to right (first exon to last exon)
						for my $k (0 .. @exonstart-1) {
							$k and $lenintron += ($exonstart[$k]-$exonend[$k-1]-1);		#calculate cumulative intron length
							$lenexon += ($exonend[$k]-$exonstart[$k]+1);
							if ($cdsstart >= $exonstart[$k]) {				#calculate CDS start accurately by considering intron length
								$rcdsstart = $cdsstart-$txstart-$lenintron+1;
								
								if ($cdsstart <= $exonend[$k]) {	#CDS start is within this exon
									$lenexon = ($exonend[$k]-$cdsstart+1);
								} else {				#CDS start is in previous exon
									#$lenexon += ($exonend[$k]-$exonstart[$k]+1);
								}
							}
							
							#splicing calculation
							if ($start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1 or $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
								$splicing{$name2}++;		#when query start site is close to exon start or exon end
							}
							if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1 or $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
								$splicing{$name2}++;		#when query end site is close to exon start or exon end
							}
							if ($start <= $exonstart[$k] and $end>=$exonstart[$k] or $start <= $exonend[$k] and $end >= $exonend[$k]) {
								$splicing{$name2}++;		#when query encompass the exon/intron boundary
							}
							
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if this splicing leads to amino acid change (rather than UTR change)
							if ($splicing{$name2} and $start==$end and $start>=$cdsstart) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$lenexon -= ($exonend[$k]-$exonstart[$k]);		#formerly "$exonend[$k]-$exonstart[$k]+1"; changed this line 2011oct01 given German's comments. The coding portion for acceptor site should be added by one, but for donor site should be fine
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:c.$lenexon-" . ($exonstart[$k]-$start) . "$ref>$obs,";
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:c.$lenexon+" . ($start-$exonend[$k]) . "$ref>$obs,";
								}
							}
							
							if ($start < $exonstart[$k]) {
								if ($end >= $exonstart[$k]) {	#exonic 
									$rvarstart = $exonstart[$k]-$txstart-$lenintron+1;
									
									for my $m ($k .. @exonstart-1) {
										$m > $k and $lenintron += ($exonstart[$m]-$exonend[$m-1]-1);
										if ($end < $exonstart[$m]) {
											#query           --------
											#gene     <--**---******---****---->
											$rvarend = $exonend[$m-1]-$txstart-$lenintron+1 + ($exonstart[$m]-$exonend[$m-1]-1);
											last;
										} elsif ($end <= $exonend[$m]) {
											#query           -----------
											#gene     <--**---******---****---->
											$rvarend = $end-$txstart-$lenintron+1;
											last;
										}
									}
									if (not defined $rvarend) {
										$rvarend = $txend-$txstart-$lenintron+1;		#if this value is longer than transcript length, it suggest whole gene deletion
									}
									
									#here the trick begins to differentiate UTR versus coding exonic
									if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
										#query  ----
										#gene     <--*---*->
										$utr5{$name2}++;		#positive strand for UTR5
									} elsif ($start > $cdsend) {
										#query             ----
										#gene     <--*---*->
										$utr3{$name2}++;		#positive strand for UTR3
									} else {									
										$exonic{$name2}++;
										not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '+', $i, $k+1, $nextline];	#refseq CDS start, refseq variant start. obs is non-zero (obs is specified by user)
									}
									$foundgenic++;
									last;
								} elsif ($k and $start > $exonend[$k-1]) {	#intronic
									$intronic{$name2}++;
									$foundgenic++;
									last;
								}
							} elsif ($start <= $exonend[$k]) {	#exonic
								$rvarstart = $start-$txstart-$lenintron+1;
								
								for my $m ($k .. @exonstart-1) {
									$m > $k and $lenintron += ($exonstart[$m]-$exonend[$m-1]-1);
									if ($end < $exonstart[$m]) {
										#query              ------
										#gene     <--**---******---****---->
										$rvarend = $exonend[$m-1]-$txstart-$lenintron+1 + ($exonstart[$m]-$exonend[$m-1]-1);
										last;
									} elsif ($end <= $exonend[$m]) {
										#query           -----------
										#gene     <--**---******---****---->
										$rvarend = $end-$txstart-$lenintron+1;
										last;
									}
								}
								if (not defined $rvarend) {
									$rvarend = $txend-$txstart-$lenintron+1;		#if this value is longer than transcript length, it suggest whole gene deletion
								}
								
								#here is the trick begins to differentiate UTR versus coding exonic
								if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
									#query  ----
									#gene     <--*---*->
									$utr5{$name2}++;		#positive strand for UTR5
								} elsif ($start > $cdsend) {
									#query             ----
									#gene     <--*---*->
									$utr3{$name2}++;		#positive strand for UTR3
								} else {
									$exonic{$name2}++;
									not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '+', $i, $k+1, $nextline];		#queryindex, refseq CDS start, refseq variant start
								}
								$foundgenic++;
								last;
							}
						}
					} elsif ($dbstrand eq '-') {		#process negative strand (in the future, this should be fused to the paragraph above for positive strands; for now, I keep them separate for easier debugging)
						for (my $k = @exonstart-1; $k>=0; $k--) {
							$k < @exonstart-1 and $lenintron += ($exonstart[$k+1]-$exonend[$k]-1);
							$lenexon += ($exonend[$k]-$exonstart[$k]+1);
							if ($cdsend <= $exonend[$k]) {		#calculate CDS start accurately by considering intron length
								$rcdsstart = $txend-$cdsend-$lenintron+1;
								
								if ($cdsend >= $exonstart[$k]) {	#CDS start within this exon
									$lenexon = ($cdsend-$exonstart[$k]+1);
								} else {				#CDS start in prevous exon
									#$lenexon += ($exonend[$k]-$exonstart[$k]+1);
								}
							}
							
							#splicing calculation
							if ($start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1 or $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
								$splicing{$name2}++;
							}
							if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1 or $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
								$splicing{$name2}++;
							}
							if ($start <= $exonstart[$k] and $end>=$exonstart[$k] or $start <= $exonend[$k] and $end >= $exonend[$k]) {
								$splicing{$name2}++;
							}
							
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if this splicing leads to amino acid change (rather than UTR change)
							if ($splicing{$name2} and $start==$end and $start<=$cdsend) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-$k+1)}:c.$lenexon+" . ($exonstart[$k]-$start) . revcom($ref) . '>' . revcom ($obs) . ',';
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$lenexon -= ($exonend[$k]-$exonstart[$k]);	#formerly "$exonend[$k]-$exonstart[$k]+1"; changed this line 2011oct01 given German's comments. The coding portion for acceptor site should be added by one, but for donor site should be fine
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-$k+1)}:c.$lenexon-" . ($start-$exonend[$k]) . revcom($ref) . '>' . revcom($obs) . ',';
								}
							}
							
							if ($end > $exonend[$k]) {
								if ($start <= $exonend[$k]) {
									$rvarstart = $txend-$exonend[$k]-$lenintron+1;
									
									for (my $m = $k; $m >= 0; $m--) {
										$m < $k and $lenintron += ($exonstart[$m+1]-$exonend[$m]-1);
										if ($start > $exonend[$m]) {
											#query           --------
											#gene     <--**---******---****---->
											#$rvarend = $txend-$exonstart[$m]-$lenintron+1 - ($exonstart[$m+1]-$exonend[$m]-1);	#commented out 2011feb18
											$rvarend = $txend-$exonstart[$m+1]+1-$lenintron + ($exonstart[$m+1]-$exonend[$m]-1);	#fixed this 2011feb18
											last;		#finsih the cycle!!!!!!!!!!!!!!!!!!!
										} elsif ($start >= $exonstart[$m]) {		#start within exons
											#query               ----
											#gene     <--**---******---****---->
											$rvarend = $txend-$start-$lenintron+1;
											last;
										}
									}
									if (not defined $rvarend) {				#if rvarend is not found, then the whole tail of gene is covered
										$rvarend = $txend-$txstart-$lenintron+1;
									}
									
									#here is the trick begins to differentiate UTR versus coding exonic
									if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
										#query  ----
										#gene     <--*---*->
										$utr3{$name2}++;		#negative strand for UTR5
									} elsif ($start > $cdsend) {
										#query             ----
										#gene     <--*---*->
										$utr5{$name2}++;		#negative strand for UTR3
									} else {
										$exonic{$name2}++;
										not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '-', $i, @exonstart-$k, $nextline];
									}
									$foundgenic++;
									last;
								} elsif ($k < @exonstart-1 and $end < $exonstart[$k+1]) {
									$intronic{$name2}++;
									$foundgenic++;
									last;
								}
							} elsif ($end >= $exonstart[$k]) {
								$rvarstart = $txend-$end-$lenintron+1;		#all the rvarstart, rvarend are with respect to the cDNA sequence (so rvarstart corresponds to end of variants)
								
								for (my $m = $k; $m >= 0; $m--) {
									$m < $k and $lenintron += ($exonstart[$m+1]-$exonend[$m]-1);
									if ($start > $exonend[$m]) {
										#query           ----
										#gene     <--**---******---****---->
										#$rvarend = $txend-$exonstart[$m]-$lenintron+1 - ($exonstart[$m+1]-$exonend[$m]-1);		#commented out 2011feb18 due to bug (10 42244567 42244600 CACCTTTGCTTGATATGATAATATAGTGCCAAGG - hetero)
										$rvarend = $txend-$exonstart[$m+1]+1 - $lenintron + ($exonstart[$m+1]-$exonend[$m]-1);		#fixed this 2011feb18
										last;			#finish the circle of counting exons!!!!!
									} elsif ($start >= $exonstart[$m]) {			#the start is right located within exon
										#query        -------
										#gene     <--**---******---****---->
										$rvarend = $txend-$start-$lenintron+1;
										last;						#finish the cycle
									}
								}
								if (not defined $rvarend) {					#if rvarend is not found, then the whole tail of gene is covered
									$rvarend = $txend-$txstart-$lenintron+1;
								}
								
								#here the trick begins to differentiate UTR versus coding exonic
								if ($end < $cdsstart) {			#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
									#query  ----
									#gene     <--*---*->
									$utr3{$name2}++;		#negative strand for UTR5
								} elsif ($start > $cdsend) {
									#query             ----
									#gene     <--*---*->
									$utr5{$name2}++;		#negative strand for UTR3
								} else {
									$exonic{$name2}++;
									not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '-', $i, @exonstart-$k, $nextline];
								}
								$foundgenic++;
								last;
							}
						}
					}
				}
			}
		}
		$foundgenic or $intergenic{''}++;		#changed $name2 to '' on 20110924
		$i =~ m/000000$/ and printerr "NOTICE: Finished analyzing $i query variants\n";

	
		my (@txname, %genename);
		my (%newsplicing);
		
		#process splicing annotation (change %splicing hash to %newsplicing, where gene name is replaced by gene name plus splicing annotation)
		if (%splicing) {
			if ($end-$start+1<=$splicing_threshold) {		#make sure that long indel are not considered here
				for my $tempname (keys %splicing) {
					if ($splicing_anno{$tempname}) {
						$splicing_anno{$tempname} =~ s/,$//;	#remove the trailing comma
						$tempname .= "($splicing_anno{$tempname})";
					}
					$newsplicing{$tempname}++;
				}
			} else {
				%newsplicing = %splicing;
			}
		}
		
		
		if ($separate) {	#separately print out each effect on one line
			if (%exonic or %splicing or %intronic or %utr5 or %utr3 or %ncrna or %upstream or %downstream) {
				#if (%ncrna) {
				#	%exonic and print OUT "ncRNA_exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#	%splicing and $end-$start+1<=$splicing_threshold and print OUT "ncRNA_splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#	%intronic and print OUT "ncRNA_intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
				#	
				#	for my $key (keys %ncrna) {
				#		delete $exonic{$key};
				#		delete $splicing{$key};
				#		delete $intronic{$key};
				#	}
				#}
				if (%exonic) {
					my (@coding, @noncoding);
					for my $key (keys %exonic) {
						if ($ncrna{$key}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
					@coding and print OUT "exonic\t", join(",", sort @coding), "\t", $nextline, "\n";
					@noncoding and print OUT "ncRNA_exonic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
				if (%splicing and $end-$start+1<=$splicing_threshold) {
					my (@coding, @noncoding);
					for my $key (keys %splicing) {
						if ($ncrna{$key}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
					@coding and print OUT "splicing\t", join(",", sort @coding), "\t", $nextline, "\n";
					@noncoding and print OUT "ncRNA_splicing\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
				if (%intronic) {
					my (@coding, @noncoding);
					# print Dumper (\%intronic);<STDIN>;
					for my $key (keys %intronic) {
						if ($ncrna{$key}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
					@coding and print OUT "intronic\t", join(",", sort @coding), "\t", $nextline, "\n";
					@noncoding and print OUT "ncRNA_intronic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
					print INTRON "$nextline\n";
				}
				
				#the following paragraph is commented out on 2011oct02
				#%exonic and print OUT "exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#%splicing and $end-$start+1<=$splicing_threshold and print OUT "splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#%intronic and print OUT "intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
				%utr5 and print OUT "UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
				%utr3 and print OUT "UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
				#if (%ncrna) {
				#	if (%exonic) {
				#		print OUT "ncRNA_exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#	}
				#	if (%splicing and $end-$start+1<=$splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
				#		print OUT "ncRNA_splicing\t",join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				#	}
				#	if (%utr5) {	#ncRNA should not have UTR5 or UTR3. If such an output exists, then there is a bug that should be reported and debugged!!!!!!!!!!!!!!!!
				#		print OUT "ncRNA_UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
				#	}
				#	if (%utr3) {	#ncRNA should not have UTR5 or UTR3. If such an output exists, then there is a bug that should be reported and debugged!!!!!!!!!!!!!!!!
				#		print OUT "ncRNA_UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
				#	}
				#	if (%intronic) {
				#		print OUT "ncRNA_intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
				#	}
				#}
				%upstream and print OUT "upstream\t", join(",", sort keys %upstream), "\t", $nextline, "\n";
				%downstream and print OUT "downstream\t", join(",", sort keys %downstream), "\t", $nextline, "\n";
				print INTRON "$nextline\n" if (%upstream || %downstream);
			} elsif (%intergenic) {
				$genel ||= "NONE";
				$gener ||= "NONE";
				$distl ||= "NONE";
				$distr ||= "NONE";
				print OUT "intergenic\t", "$genel(dist=$distl),$gener(dist=$distr)", "\t", $nextline, "\n";
				print INTRON "$nextline\n";
			} else {
				die "FATAL ERROR: please report bug to ANNOVAR author with your input file\n";
			}
		} else {			
			if (@precedence) {
				my $foundmatch;
				for my $i (0 .. @precedence-2) {
					$precedence[$i] eq 'exonic' and %exonic and $foundmatch++;
					$precedence[$i] eq 'splicing' and %splicing and $foundmatch++;
					$precedence[$i] eq 'intronic' and %intronic and $foundmatch++;
					$precedence[$i] eq 'utr5' and %utr5 and $foundmatch++;
					$precedence[$i] eq 'utr3' and %utr3 and $foundmatch++;
					$precedence[$i] eq 'ncrna' and %ncrna and $foundmatch++;
					$precedence[$i] eq 'upstream' and %upstream and $foundmatch++;
					$precedence[$i] eq 'downstream' and %downstream and $foundmatch++;
					$precedence[$i] eq 'intergenic' and %intergenic and $foundmatch++;
					if ($foundmatch) {
						for my $j ($i+1 .. @precedence-1) {
							$precedence[$j] eq 'exonic' and %exonic = ();
							$precedence[$j] eq 'splicing' and %splicing = ();
							$precedence[$j] eq 'intronic' and %intronic = ();
							$precedence[$j] eq 'utr5' and %utr5 = ();
							$precedence[$j] eq 'utr3' and %utr3 = ();
							$precedence[$j] eq 'ncrna' and %ncrna = ();
							$precedence[$j] eq 'upstream' and %upstream = ();
							$precedence[$j] eq 'downstream' and %downstream = ();
							$precedence[$j] eq 'intergenic' and %intergenic = ();
						}
						last;
					}
				}
			}
			
		
				
			if (%exonic) {
				my (@coding, @noncoding);
				for my $key (keys %exonic) {
					if ($ncrna{$key}) {
						push @noncoding, $key;
					} else {
						push @coding, $key;
					}
				}
				if (@coding and %splicing and $end-$start+1<=$splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
					print OUT "exonic;splicing\t", join(",", sort @coding), ";", join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				} elsif (@coding) {
					print OUT "exonic\t", join(",", sort @coding), "\t", $nextline, "\n";
				} elsif (@noncoding) {
					print OUT "ncRNA_exonic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
			} elsif (%splicing) {
				my (@coding, @noncoding);
				for my $key (keys %newsplicing) {
					$key =~ m/^([^\(]+)/;
					if ($ncrna{$1}) {
						push @noncoding, $key;
					} else {
						push @coding, $key;
					}
				}
				if (@coding) {
					print OUT "splicing\t", join (',', sort @coding), "\t", $nextline, "\n";
				} elsif (@noncoding) {
					print OUT "ncRNA_splicing\t", join (',', sort @noncoding), "\t", $nextline, "\n";
				}
			} elsif (%ncrna) {
				my (@coding, @noncoding);
				for my $key (keys %intronic) {
					if ($ncrna{$key}) {
						push @noncoding, $key;
					} else {
						push @coding, $key;
					}
				}
				#print OUT "ncRNA\t", join(",", sort keys %ncrna), "\t", $nextline, "\n";
				
				#if (%exonic) {
				#	if (%splicing and $end-$start+1<=$splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
				#		print OUT "ncRNA_exonic;ncRNA_splicing\t", join(",", sort keys %exonic), ";", join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				#	} else {
				#		print OUT "ncRNA_exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#	}
				#} elsif (%splicing) {
				#	print OUT "ncRNA_splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#} elsif (%utr5 or %utr3) {		#ncRNA should not have UTR5 or UTR3. If such an output exists, then there is a bug that should be reported and debugged!!!!!!!!!!!!!!!!
				if (%utr5 or %utr3) {
					if (%utr5 and %utr3) {
						print OUT "ncRNA_UTR5;ncRNA_UTR3\t", join(",", sort keys %utr5), ";", join(",", sort keys %utr3), "\t", $nextline, "\n";		#use ";" to separate UTR5 and UTR3 genes
					} elsif (%utr5) {
						print OUT "ncRNA_UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
					} else {
						print OUT "ncRNA_UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
					}
				} elsif (@noncoding) {
					print OUT "ncRNA_intronic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				} else {
					die "FATAL ERROR: please report bug to ANNOVAR author with your input line <$nextline>\n";
				}
			} elsif (%utr5 or %utr3) {
				if (%utr5 and %utr3) {
					print OUT "UTR5;UTR3\t", join(",", sort keys %utr5), ";", join(",", sort keys %utr3), "\t", $nextline, "\n";		#use ";" to separate UTR5 and UTR3 genes
				} elsif (%utr5) {
					print OUT "UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
				} else {
					print OUT "UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
				}
			} elsif (%intronic) {
				print OUT "intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
			} elsif (%upstream or %downstream) {
				if (%upstream and %downstream) {
					print OUT "upstream;downstream\t", join(",", sort keys %upstream), ";", join(",", sort keys %downstream), "\t", $nextline, "\n";
				} elsif (%upstream) {
					print OUT "upstream\t", join(",", sort keys %upstream), "\t", $nextline, "\n";
				} else {
					print OUT "downstream\t", join(",", sort keys %downstream), "\t", $nextline, "\n";
				}
			} elsif (%intergenic) {
				$genel ||= "NONE";
				$gener ||= "NONE";
				$distl ||= "NONE";
				$distr ||= "NONE";
				print OUT "intergenic\t", "$genel(dist=$distl),$gener(dist=$distr)", "\t", $nextline, "\n";
			} else {
				die "FATAL ERROR: please report bug to ANNOVAR author with your input line <$nextline>\n";
			}
		} 
		print INTRON "$nextline\n" if (%intronic || %upstream || %downstream || %intergenic);
	}
	%refseqvar and annotateExonicVariants (\%refseqvar, $geneidmap, $cdslen, $mrnalen);
	return ($linecount, $invalidcount);
}

sub annotateExonicVariants {
	my ($refseqvar, $geneidmap, $cdslen, $mrnalen) = @_;
	my $refseqhash;
	my $function = {};
	my %varinfo;				#variants information (same as input line)
	my %unmatch_wtnt_ref;			#count how many user input variant has wildtype nucleotide that does not match reference genome (use hash, because if we use count, we are merely counting variant-transcript combinations)
	my @unmatch_example;
	
	$refseqhash = readSeqFromFASTADB ($refseqvar);
	for my $seqid (keys %$refseqvar) {
		for my $i (0 .. @{$refseqvar->{$seqid}}-1) {
			my ($refcdsstart, $refvarstart, $refvarend, $refstrand, $index, $exonpos, $nextline) = @{$refseqvar->{$seqid}->[$i]};
			my ($wtnt3, $wtnt3_after, @wtnt3, $varnt3, $wtaa, $wtaa_after, $varaa, $varpos);		#wtaa_after is the aa after the wtaa
			my ($chr, $start, $end, $ref, $obs);
			my @nextline = split (/\s+/, $nextline);
			($chr, $start, $end, $ref, $obs) = @nextline[@avcolumn];
			($ref, $obs) = (uc $ref, uc $obs);
			$zerostart and $start++;
			$chr =~ s/^chr//i;
			
			$varinfo{$index} = $nextline;			
			if (not $refseqhash->{$seqid}) {					#this refseq do not have FASTA sequence so cannot be interrogated
				$function->{$index}{unknown} = "UNKNOWN";
				next;
			}
			my $ref_tag='';
			my $fs = (($refvarstart-$refcdsstart) % 3);#hmv this is the phase
			if ($refvarstart-$fs-1 > length($refseqhash->{$seqid})) {
				printerr "WARNING: Potential database annotation error seqid=$seqid, refvarstart=$refvarstart, fs=$fs, seqlength=", length($refseqhash->{$seqid}), " refcdsstart=$refcdsstart, with inputline=$nextline\n";
				# next;#hv taken out 9/22/14
			}
			$wtnt3 = substr ($refseqhash->{$seqid}, $refvarstart-$fs-1, 3);
			### added hv ####
			if ($fs==0){
				$wtnt3=uc(substr($wtnt3,0,1)).lc(substr($wtnt3,1,2));
			}elsif ($fs==1){
				$wtnt3=lc(substr($wtnt3,0,1)).uc(substr($wtnt3,1,1)).lc(substr($wtnt3,2,1));
			}elsif ($fs==2){
				$wtnt3=lc(substr($wtnt3,0,2)).uc(substr($wtnt3,2,1));
			}
			
			####end hv block
			if (length ($refseqhash->{$seqid}) >= $refvarstart-$fs+3) {	#going into UTR
				$wtnt3_after = substr ($refseqhash->{$seqid}, $refvarstart-$fs+2, 3);
			} else {
				$wtnt3_after = '';					#last amino acid in the sequence without UTR (extremely rare situation) (example: 17        53588444        53588444        -       T       414     hetero)
			}
			@wtnt3 = split (//, $wtnt3);
			my $addon=":$wtnt3";
			if (@wtnt3 != 3 and $refvarstart-$fs-1>=0) {			#some times there are database annotation errors (example: chr17:3,141,674-3,141,683), so the last coding frame is not complete and as a result, the cDNA sequence is not complete
				$function->{$index}{unknown} = "UNKNOWN";
				next;
			}
			
			if ($refstrand eq '-') {					#change the observed nucleotide to the reverse strand
				$obs = revcom ($obs);
				$ref = revcom ($ref);
			}
			
			if ($start == $end) {
				if ($ref eq '-') {					#insertion variant
					#the insertion coordinate system in ANNOVAR always uses "position after the current site"
					#in positive strand, this is okay
					#in negative strand, the "after current site" becomes "before current site" during transcription
					#therefore, appropriate handling is necessary to take this into account
					#for example, for a trinucleotide GCC with frameshift of 1 and insertion of CCT
					#in positive strand, it is G-CTT-CC
					#but if the transcript is in negative strand, the genomic sequence should be GC-CCT-C, and transcript is G-AGG-GC
					if ($refstrand eq '+') {
						if ($fs == 1) {
							$varnt3 = $wtnt3[0] . $wtnt3[1] . $obs . $wtnt3[2];
						} elsif ($fs == 2) {
							$varnt3 = $wtnt3[0] . $wtnt3[1] . $wtnt3[2] . $obs;
						} else {
							$varnt3 = $wtnt3[0] . $obs . $wtnt3[1] . $wtnt3[2];
						}
					} elsif ($refstrand eq '-') {
						if ($fs == 1) {
							$varnt3 = $wtnt3[0] . $obs . $wtnt3[1] . $wtnt3[2];
						} elsif ($fs == 2) {
							$varnt3 = $wtnt3[0] . $wtnt3[1] . $obs . $wtnt3[2];
						} else {
							$varnt3 = $obs . $wtnt3[0] . $wtnt3[1] . $wtnt3[2];
						}
					}
					#hv note: foreach function, added $wtnt3 to the end of the impact assessment
					($wtaa, $wtaa_after, $varaa, $varpos) = (translateDNA ($wtnt3), translateDNA ($wtnt3_after), translateDNA ($varnt3), int(($refvarstart-$refcdsstart)/3)+1);
					$wtaa_after and $wtaa_after eq '*' and $wtaa_after = 'X';		#wtaa_after could be undefined, if the current aa is the stop codon (X) (example: 17        53588444        53588444        -       T)

					my $canno = "c." . ($refvarstart-$refcdsstart+1) .  "_" . ($refvarstart-$refcdsstart+2) . "ins$obs";		#cDNA level annotation
					if (! $relpos){ $addon='';}
					if (length ($obs) % 3 == 0) {
						if ($wtaa eq '*') {			#mutation on stop codon
							if ($varaa =~ m/\*/) {
								$varaa =~ s/\*.*/X/;	#delete all aa after stop codon, but keep the aa before
								$function->{$index}{nfsins} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "delins$varaa$addon,";		#stop codon is stil present
							} else {
								$function->{$index}{stoploss} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "delins$varaa$addon,";	#stop codon is lost
							}
						} else {
							if ($varaa =~ m/\*/) {
								$varaa =~ s/\*.*/X/;	#delete all aa after stop codon, but keep the aa before
								$function->{$index}{stopgain} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "delins$varaa$addon,";
							} else {
								$function->{$index}{nfsins} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "delins$varaa$addon,";
							}
						}
					} else {
						if ($wtaa eq '*') {			#mutation on stop codon
							if ($varaa =~ m/\*/) {		#in reality, this cannot be differentiated from non-frameshift insertion, but we'll still call it frameshift
								$varaa =~ s/\*.*/X/;	#delete all aa after stop codon, but keep the aa before
								$function->{$index}{fsins} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "delins$varaa$addon,";
							} else {
								$function->{$index}{stoploss} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "delins$varaa$addon,";
							}
						} else {
							if ($varaa =~ m/\*/) {
								$varaa =~ s/\*.*/X/;	#delete all aa after stop codon, but keep the aa before
								$function->{$index}{stopgain} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "_$wtaa_after" . ($varpos+1) . "delins$varaa$addon,";
							} else {
								$function->{$index}{fsins} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "fs$addon,";
							}
						}
					}
				} elsif ($obs eq '-') {					#single nucleotide deletion
					my $deletent;
					if ($fs == 1) {
						$deletent = $wtnt3[1];
						$varnt3 = $wtnt3[0].$wtnt3[2].$wtnt3_after;
					} elsif ($fs == 2) {
						$deletent = $wtnt3[2];
						$varnt3 = $wtnt3[0].$wtnt3[1].$wtnt3_after;
					} else {
						$deletent = $wtnt3[0];
						$varnt3 = $wtnt3[1].$wtnt3[2].$wtnt3_after;
					}
					($wtaa, $varaa, $varpos) = (translateDNA ($wtnt3), translateDNA ($varnt3),  int(($refvarstart-$refcdsstart)/3)+1);
					my $canno = "c." . ($refvarstart-$refcdsstart+1) . "del$deletent";
					if ($wtaa eq '*') {				#mutation on stop codon
						if ($varaa =~ m/\*/) {			#stop codon is still stop codon
							$function->{$index}{nfsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "X$addon,";	#changed fsdel to nfsdel on 2011feb19
						} else {				#stop codon is lost
							$function->{$index}{stoploss} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "$varaa$addon,";
						}
					} else {
						if ($varaa =~ m/\*/) {			#new stop codon created
							$function->{$index}{stopgain} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "X$addon,";
						} else {
							$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "fs$addon,";
						}
					}
				} elsif (length ($obs) > 1) {				#block substitution (since start==end, this changed from 1nt to several nt)
					if (($refvarend-$refvarstart+1-length($obs)) % 3 == 0) {
						$function->{$index}{nfssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "delins$obs$addon,";
					} else {
						$function->{$index}{fssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "delins$obs$addon,";
					}
				} else {						#single nucleotide substitution variant
					my ($canno);
					#hv added next block to accomodate possible incorrect data
					# if ($correctAllele){
					# 	if ($ref and $wtnt3[$fs] ne $ref) {
					# 		my $rev_ref=revcom($ref);
					# 		if ($rev_ref eq $wtnt3[$fs] && $refstrand eq '-'){#if the user specified the alleles with respect to the gene
					# 			$ref_tag="RC $ref->$obs(+/+)";
					# 			$ref=$rev_ref;$obs=revcom($obs);
					# 			$verbose and printerr "WARNING: ALLELE MISMATCH(revcmp): strand=$refstrand user-specified-allele=$ref annovar-inferred-allele=$wtnt3[1] inputline=<$nextline>\n";
								
					# 		}elsif ($obs  eq "$ref"){#if the user mixed up the $ref with the $obs alleles
					# 			my $tmp=$ref;$ref=$obs;$obs=$tmp;
					# 			$verbose and printerr "WARNING: ALLELE MISMATCH(ref vs mt): strand=$refstrand user-specified-allele=$ref annovar-inferred-allele=$wtnt3[1] inputline=<$nextline>\n";
					# 			$ref_tag="SW $obs->$tmp";
					# 		}else{#use default annovar behavior by ; can't determine what the user did, so ignore the reference and use $obs as is.
								$unmatch_wtnt_ref{$index}++;
								@unmatch_example or push @unmatch_example, $nextline;
								$verbose and printerr "WARNING: ALLELE MISMATCH: strand=$refstrand user-specified-allele=$ref annovar-inferred-allele=$wtnt3[1] inputline=<$nextline>\n";
								# $ref_tag="IR $ref and $obs,$wtnt3[$fs] ne $ref";
					# 		}
					# 	}
					# }
					#end hv
					if ($fs == 1) {
						$varnt3 = $wtnt3[0] . $obs . $wtnt3[2];
						$canno = "c.$wtnt3[1]" . ($refvarstart-$refcdsstart+1) . $obs;
						$hgvs and $canno = "c." . ($refvarstart-$refcdsstart+1) . $wtnt3[1] . ">" . $obs;
						if ($ref and $wtnt3[1] ne $ref) {
							$unmatch_wtnt_ref{$index}++;
							@unmatch_example or push @unmatch_example, $nextline;
							$verbose and printerr "WARNING: ALLELE MISMATCH: strand=$refstrand user-specified-allele=$ref annovar-inferred-allele=$wtnt3[1] inputline=<$nextline>\n";
						}
					} elsif ($fs == 2) {
						$varnt3 = $wtnt3[0] . $wtnt3[1]. $obs;
						$canno = "c.$wtnt3[2]" . ($refvarstart-$refcdsstart+1) . $obs;
						$hgvs and $canno = "c." . ($refvarstart-$refcdsstart+1) . $wtnt3[2] . ">" . $obs;
						if ($ref and $wtnt3[2] ne $ref) {
							$unmatch_wtnt_ref{$index}++;
							@unmatch_example or push @unmatch_example, $nextline;
							$verbose and printerr "WARNING: ALLELE MISMATCH: strand=$refstrand user-specified-allele=$ref annovar-inferred-allele=$wtnt3[2] inputline=<$nextline>\n";
						}
					} else {
						$varnt3 = $obs . $wtnt3[1] . $wtnt3[2];
						$canno = "c.$wtnt3[0]" . ($refvarstart-$refcdsstart+1) . $obs;
						$hgvs and $canno = "c." . ($refvarstart-$refcdsstart+1) . $wtnt3[0] . ">" . $obs;
						if ($ref and $wtnt3[0] ne $ref) {
							$unmatch_wtnt_ref{$index}++;
							@unmatch_example or push @unmatch_example, $nextline;
							$verbose and printerr "WARNING: ALLELE MISMATCH: strand=$refstrand user-specified-allele=$ref annovar-inferred-allele=$wtnt3[0] inputline=<$nextline>\n";
						}
					}
					($wtaa, $varaa, $varpos) = (translateDNA ($wtnt3), translateDNA ($varnt3), int(($refvarstart-$refcdsstart)/3)+1);
					if (!$correctAllele){
						$ref_tag='';
					}else{$ref_tag.="\t";}
					if ($wtaa eq $varaa) {
						$wtaa eq '*' and ($wtaa, $varaa) = qw/X X/;		#change * to X in the output
						$function->{$index}{ssnv} .= "$ref_tag$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos$varaa$addon,";
					} elsif ($varaa eq '*') {
						$function->{$index}{stopgain} .= "$ref_tag$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa${varpos}X$addon,";
					} elsif ($wtaa eq '*') {
						$function->{$index}{stoploss} .= "$ref_tag$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos$varaa$addon,";
					} else {
						$function->{$index}{nssnv} .= "$ref_tag$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos$varaa$addon,";
					}
				}
			} elsif ($obs eq '-') {				#deletion variant involving several nucleotides
				($wtaa, $varpos) = (translateDNA ($wtnt3), int(($refvarstart-$refcdsstart)/3)+1);		#wildtype amino acid, position of amino acid
				my ($varposend, $canno);		#the position of the last amino acid in the deletion
				if ($refvarstart<=$refcdsstart) {	#since the first amino acid is deleted, the whole gene is considered deleted
					$function->{$index}{stdel} .= "$geneidmap->{$seqid}:$seqid:wholegene$addon,";	#it is exonic variant, so the varend has to hit the first exon
				} elsif ($refvarend >= $cdslen->{$seqid}+$refcdsstart) {	#3' portion of the gene is deleted
					$varposend = int ($cdslen->{$seqid}/3);		#cdslen should be multiples of 3, but just in case of database mis-annotation
					$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($cdslen->{$seqid}+$refcdsstart-1) . "del";
					$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.${varpos}_${varposend}del$addon,";
				} elsif (($refvarend-$refvarstart+1) % 3 == 0) {
					$varposend = int (($refvarend-$refcdsstart)/3) + 1;
					$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "del";
					$function->{$index}{nfsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.${varpos}_${varposend}del$addon,";
				} else {
					$varposend = int (($refvarend-$refcdsstart)/3) + 1;
					$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "del";
					$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.${varpos}_${varposend}del$addon,";
				}
			} else {							#block substitution event
				if (($refvarend-$refvarstart+1-length($obs)) % 3 == 0) {
					$function->{$index}{nfssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "$obs$addon,";
				} else {
					$function->{$index}{fssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "$obs$addon,";
				}
			}
		}
	}
	for my $index (sort {$a<=>$b} keys %$function) {
		if ($separate) {			#print out each type of exonic mutations separately (one effect in one line), rather than printing out only the most important function
			if ($function->{$index}{fsins}) {
				print EXONIC "line$index\t", "frameshift insertion\t$function->{$index}{fsins}\t", $varinfo{$index}, "\n";
			}
			if ($function->{$index}{fsdel}) {
				print EXONIC "line$index\t", "frameshift deletion\t$function->{$index}{fsdel}\t", $varinfo{$index}, "\n";
			}
			if ($function->{$index}{stdel}) {
				print EXONIC "line$index\t", "start codon deletion\t$function->{$index}{fsdel}\t", $varinfo{$index}, "\n";
			}
			if ($function->{$index}{fssub}) {
				print EXONIC "line$index\t", "frameshift substitution\t$function->{$index}{fssub}\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{stopgain}) {
				print EXONIC "line$index\t", "stopgain SNV\t$function->{$index}{stopgain}\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{stoploss}) {
				print EXONIC "line$index\t", "stoploss SNV\t$function->{$index}{stoploss}\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{nfsins}) {
				print EXONIC "line$index\t", "nonframeshift insertion\t$function->{$index}{nfsins}\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{nfsdel}) {
				print EXONIC "line$index\t", "nonframeshift deletion\t$function->{$index}{nfsdel}\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{nfssub}) {
				print EXONIC "line$index\t", "nonframeshift substitution\t$function->{$index}{nfssub}\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{nssnv}) {
				print EXONIC "line$index\t", "nonsynonymous SNV\t$function->{$index}{nssnv}\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{ssnv}) {
				print EXONIC "line$index\t", "synonymous SNV\t$function->{$index}{ssnv}\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{unknown}) {
				print EXONIC "line$index\t", "unknown\t$function->{$index}{unknown}\t", $varinfo{$index}, "\n";
			}
		} else {				#print out only the most important functional changes (for example, chr3:9931279-9931279 G->A can be both non-synonymous and synonymous mutations based on UCSC gene model)
			print EXONIC "line$index\t";
			my $sortout;
			if ($sortout = $function->{$index}{fsins}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "frameshift insertion\t$sortout\t";
			} elsif ($sortout = $function->{$index}{fsdel}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "frameshift deletion\t$sortout\t";
			} elsif ($sortout = $function->{$index}{fssub}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "frameshift substitution\t$sortout\t";
			} elsif ($sortout = $function->{$index}{stopgain}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "stopgain SNV\t$sortout\t";
			} elsif ($sortout = $function->{$index}{stoploss}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "stoploss SNV\t$sortout\t";
			} elsif ($sortout = $function->{$index}{nfsins}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "nonframeshift insertion\t$sortout\t";
			} elsif ($sortout = $function->{$index}{nfsdel}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "nonframeshift deletion\t$sortout\t";
			} elsif ($sortout = $function->{$index}{nfssub}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "nonframeshift substitution\t$sortout\t";
			} elsif ($sortout = $function->{$index}{nssnv}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "nonsynonymous SNV\t$sortout\t";
			} elsif ($sortout = $function->{$index}{ssnv}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "synonymous SNV\t$sortout\t";
			} elsif ($sortout = $function->{$index}{unknown}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "unknown\t$sortout\t";
			} elsif ($sortout = $function->{$index}{stdel}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				print EXONIC "start codon deletion\t$sortout\t";
			}
			print EXONIC $varinfo{$index}, "\n";
		}
	}
	
	if (%unmatch_wtnt_ref) {
		my $show=($#unmatch_example<10)?$#unmatch_example:10;
		printerr "\n$show?\n";
		printerr "---------------------------------------------------------------------------------------\n";
		printerr "WARNING: ${\(scalar keys %unmatch_wtnt_ref)} exonic SNPs have WRONG reference alleles specified in your input file!\n";
		printerr "WARNING: An example input line is <". join ("\n",@unmatch_example[0-$show]). ">\n";
		printerr "WARNING: ANNOVAR can still annotate exonic_variant_function for the mutation correctly!\n";
		printerr "WARNING: you may have used wrong -buildver, or did not specify the correct reference allele!\n";
		printerr "---------------------------------------------------------------------------------------\n\n";
	}
}

sub sortExonicAnnotation {
	my ($anno) = @_;
	my @anno1 = split (/,/, $anno);
	my @anno2;
	#example: ATG16L1:NM_198890:exon5:c.A409G:p.T137A,ATG16L1:NM_017974:exon8:c.A841G:p.T281A
	for my $i (0 .. @anno1-1) {
		my @temp = split (/:/, $anno1[$i]);
		if ($temp[2] =~ s/^exon//) {		#some are wholegene, not exon
			push @anno2, [$anno1[$i], @temp];
		}
	}
	if (@anno2 and @anno1==@anno2) {
		@anno2 = sort {$a->[3] <=> $b->[3] or $a->[2] cmp $b->[2]} @anno2;		#first sort by exon number, then by transcript name
		my @anno3 = map {$_->[0]} @anno2;
		return join (',', @anno3);
	} else {
		return $anno;
	}
}

sub filterQuery {
	my $xname=$dbtype1;
	if ($dbtype1=~/generic/){
		$xname=$genericdbfile;
		$xname=~s/(^$buildver\_|.txt)//g;
		# $xname.=".exact";
	}
	open (FIL, ">$outfile.${buildver}_${xname}_filtered") or die "Error: cannot write to output file $outfile.${buildver}_${xname}_filtered: $!\n"; 
	open (DROPPED, ">$outfile.${buildver}_${xname}_dropped") or die "Error: cannot write to output file $outfile.${buildver}_${xname}_dropped: $!\n";
	if ($keepline && !$noheader && !$headersWanted){
		print FIL "0LineNbr\t#$xname\n";
	}elsif($keepline && !$noheader && $headersWanted){
		my @headersWanted = split(",",$headersWanted);
		print FIL "0LineNbr\t#" . join (";",@headersWanted) . "\n";
	}
	open (INVALID, ">$outfile.invalid_input") or die "Error: cannot write to output file $outfile.invalid_input: $!\n";
	printerr "NOTICE: Variants matching filtering criteria are written to $outfile.${buildver}_${xname}_dropped, other variants are written to $outfile.${buildver}_${xname}_filtered\n";
	
	open (QUERY, $queryfile) or die "Error: cannot read from query file $queryfile: $!\n";
	
	my (%variant, $filedone, $batchdone);
	my ($linecount, $batchlinecount, $invalid, $invalidcount) = (0, 0);
	my ($chr, $start, $end, $ref, $obs, $info);
	while (1) {
		$_ = <QUERY>;
		if (not defined $_) {
			$filedone++;
		} else {
			s/[\r\n]//g;
			
			if (m/^#/ and $comment) {				#comment line start with #, do not include this is $linecount
				print FIL "#\t" if ($keepline);
				print FIL "$_\n";
				next;
			}elsif(m/^#/){
				next;
			}
			$linecount++;
			
			$batchlinecount++;
			if ($batchlinecount == $batchsize) {
				$batchdone++;
			}
			
			if ($memfree or $memtotal) {		#if these arguments are specified
				if ($linecount =~ m/00000$/) {						#about 40Mb memory per 10k lines for a typical input dataset
					my ($availmem, $allmem) = currentAvailMemory();
					$verbose and printerr "NOTICE: Current available system memory is $availmem kb (this program uses $allmem bytes memory), after reading $linecount query\n";
					if ($availmem and $availmem <= $memfree+50_000) {		#some subsequent steps may take ~50Mb memory, so here we try to allocate some more memory
						$batchdone++;
					}
					if ($memtotal and $allmem >= $memtotal-50_000) {	#when --memtotal is specified, ensure that program use less memory
						$batchdone++;
					}
				}
			}
	
			$invalid = 0;						#reset invalid status
			
			my @nextline = split (/\s+/, $_);
			($chr, $start, $end, $ref, $obs) = @nextline[@avcolumn];
			if ( not (defined $chr and defined $start and defined $end and defined $ref and defined $obs)) {
				$invalid++;
			} else {
				($ref, $obs) = (uc $ref, uc $obs);
				$zerostart and $start++;
				$chr =~ s/^chr//i;
				if ($chr =~ m/[^\w]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {
					$invalid++;
					
				} elsif ($ref eq '-' and $obs eq '-' 		#both are empty allele
					or $ref =~ m/[^ACTG0\-]/ 		#non-standard nucleotide code
					or $obs =~ m/[^ACGT0\-]/ 		#non-standard nucleotide code
					or $start =~ m/[^\d]/ 			#start is not a number
					or $end =~ m/[^\d]/ 			#end is not a number
					or $start > $end			#start is more than end
					# or $ref ne '0' and $end-$start+1 != length ($ref) 	#length mismatch with ref
					or $ref eq '-' and $start != $end	#length mismatch for insertion
					) {		
					$invalid++;

				}
			}				
			
			if ($invalid) {
				print INVALID $_, "\n";	#invalid record found
				print FIL "$linecount\tERR\n" if ($keepline);
				$invalidcount++;
				next;
			}
			if ($start == $end and $ref eq '-') {	#insertion
				$obs = "0$obs";
			} elsif ($obs eq '-') {			#deletion
				# $obs = $end-$start+1;##hvfixed 20160225 for avsnp
			} elsif ($end>$start or $start==$end and length($obs)>1) {	#block substitution	#fixed the bug here 2011feb19
				# $obs = ($end-$start+1) . $obs;
			}
			if (exists $variant{$chr, $start, $obs}) {
				$variant{$chr, $start, $obs} .= "\n$linecount\t$_";
			} else {
				$variant{$chr, $start, $obs} = "$ref\n$linecount\t$_";
			}
			
		}
		if ($filedone or $batchdone) {
			printerr "NOTICE: Processing next batch with ${\(scalar keys %variant)} unique variants in $batchlinecount input lines\n";
			filterNextBatch (\%variant);
			%variant = ();
			$batchlinecount = 0;				#reset the line count for this batch
			$batchdone = 0;
		}
		if ($filedone) {
			last;
		}
	}
	close (INVALID); close (DROPPED); close (FIL);
	if ($invalidcount) {
		printerr "NOTICE: Variants with invalid input format were written to $outfile.invalid_input\n";
	} else {
		unlink ("$outfile.invalid_input");
	}
	return ("$outfile.${buildver}_${xname}");
}

sub filterNextBatch {
	my ($variant) = @_;
	my $dbfile;
	if ($dbtype1 eq 'generic') {
		$dbfile = File::Spec->catfile ($dbloc, $genericdbfile);
	} elsif ($dbtype1 eq 'vcf') {
		$dbfile = File::Spec->catfile ($dbloc, $vcfdbfile);
	} else {
		print "Reading $dbloc/${buildver}_$dbtype1.txt\n";
		$dbfile = File::Spec->catfile ($dbloc, "${buildver}_$dbtype1.txt");
	}
	my (@record, $chr, $start, $end, $ref, $obs, $score, $qual, $fil, $info);
	my ($rsid, $strand, $ucscallele, $twoallele, $class, $af, $attribute);
	my $count_invalid_dbline;

	my ($BIN, $DBSIZE) = (0, 0);
	my %index = ();
	my $bb = {};			#a subset of %index, which corresponds to the input variants
	my $flag_idx_search = 0;	#indicate if index-based search algorithm is used (faster speed for a small number of input variants
	
	if ( -f "$dbfile.idx" ) {
		open(IDX, "$dbfile.idx") or die "Error: cannot read from input database index $dbfile.idx: $!\n";
		print 'Reading Index ' . "$dbfile.idx\n";
		my $line = <IDX>;
		if (not $line =~ m/BIN\t(\d+)\t(\d+)/) {
			printerr "WARNING: Malformed database index file $dbfile.idx.\n";
		} elsif ($2 != -s $dbfile) {		#file version is different, do not use index file in this case
			printerr "WARNING: Your index file $dbfile.idx is out of date and will not be used. ANNOVAR can still generate correct results without index file.\n";
		} else {
			($BIN, $DBSIZE) = ($1, $2);
			while ( $line = <IDX> ) {
				$line =~ s/[\r\n]+$//;
				my ( $chrom, $pos, $offset0, $offset1 ) = split (/\t/, $line);
				$chrom=~s/^chr//i;
				defined $offset1 or next;				#invalid input line in the index file
				my $quart=int(($offset1-$offset0)/4);
				$offset0=$quart if ($offset0==0);# this is for the first offset =0; do not want a negative number
				$index{"$chrom\t$pos"} = [($offset0-$quart), $offset1];
			}
		}
		close (IDX);
		print "Done with index\n";
		if (not %index) {
			printerr "WARNING: Unable to load database index successfully from file $dbfile.idx.\n";
		} else {
			foreach my $k ( keys %$variant ) {
				#my ( $chrom, $pos ) = split /\t/, $outer->{$k};
				#$chrom =~ s/.+\n//;
				my ($chrom, $pos) = split ($;, $k);
				$chrom=~s/^chr//i;
				my $bin = $pos - ( $pos % $BIN );
				if (!defined $index{"$chrom\t$bin"} ){
					my $fail=0; #keep track of failures so I don't keep looping
					until (defined  $index{"$chrom\t$bin"} || $bin<0 || $fail >10 ){
						$fail++;
						$bin-=$BIN;
					}
				}
				$bb->{ "$chrom\t$bin" }{ $k } = $variant->{ $k };
			}
			print "Done with mapping variants to indexed bins\n";
			if (scalar (keys %$bb) / scalar (keys %index) <= $indexfilter_threshold) {##HV201305
				$flag_idx_search++;
				printerr "NOTICE: Database index loaded. Total number of bins is ".  scalar (keys %index) . " and " . " number of bins to be scanned is " . scalar (keys %$bb) . "\n";
			}
		}
	} 
	if (! $flag_idx_search) {
		$bb = {1, [0, -s "$dbfile"]};
		%index = (1, [0, -s "$dbfile"]);
	}

	# print Dumper $bb;exit;
	open (DB, $dbfile) or die "Error: cannot read from input database file $dbfile: $!\n";
	printerr "NOTICE: Scanning filter database $dbfile at line#".__LINE__."\n";
# print Dumper (\%{$bb});
	my @arr;
	foreach my $b ( sort keys %$bb ) {
		next if !exists $index{$b};
		my ( $chunk_min, $chunk_max ) = @{ $index{$b} };
		#close(DB);
		#open(DB, $dbfile);
		seek( DB, $chunk_min, 0 );
		#$variant = $bb->{ $b };
		my $chunk_here = $chunk_min;
		while (<DB>) {
			my $line_length = length($_);			#calculate line length of the DB
			my (@obs2, @score2);				#for 1000G2010 data set in VCF format, some tri-allelic SNPs are present; in the future, some quad-allelic SNPs may be also present in VCF files
			s/[\r\n]+$//;
			m/\S/ or next;					#skip empty lines in the database file (sometimes this occurs)
			m/^#/ and next;					#skip the comment line
			if ($dbtype eq 'avsift') {
				@record = split (/\t/, $_);
				@record == 8 or next;#die "Error: invalid record found in DB file $dbfile (8 tab-delimited fields expected): <$_>\n";
				($chr, $start, $end, $ref, $obs, $score) = @record;
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				if ($score < $sift_threshold  && not $keepline) {		#this is a deleterious mutation, skip it (equal sign should not be used, otherwise the score=0 will be skipped)
					next;
				}
			}elsif ($dbtype=~/avsnp1\d+/) {
				@record = split (/\t/, $_);
				@record == 6 or next;#die "Error: invalid record found in DB file $dbfile (8 tab-delimited fields expected): <$_>\n";
				($chr, $start, $end, $ref, $obs,$score) = @record;
				# $score=$obs;#20151214fixed
			} elsif ($dbtype =~ m/^ljb2_/) {
				@record = split (/\t/, $_);my $pred;
				@record >= 5 or next;# die "Error: invalid record found in DB file $dbfile (at least 5 tab-delimited fields expected): <$_>\n";
				($chr, $start, $end, $ref, $obs, $score, $pred) = @record;$pred||='';
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				if (defined $score and defined $score_threshold and not $keepline) {
					if ($reverse) {
						$score > $score_threshold and next;
					} else {
						$score < $score_threshold and next;
					}
				}
				if ($dbtype=~/pp2/ && $keepline){
					if ($pred eq 'D'){
						$score="DAMAGING:$score";
					}elsif ($pred eq 'P'){
						$score="Probably DAMAGING:$score";
					}elsif ($pred eq 'B'){
						$score="Benign:$score";
					}elsif ($pred eq 'NA'){
						$score="NA:$score";
					}
				}elsif ($dbtype=~/mt/ && $keepline){
					# MutationTaster prediction, "A" ("disease_causing_automatic"), "D" ("disease_causing"), "N" ("polymorphism") or "P" ("polymorphism_automatic")
					if ($pred eq 'A'){
						$score="Disease Causing Automatic:$score";
					}elsif ($pred eq 'D'){
						$score="Disease Causing:$score";
					}elsif ($pred eq 'N'){
						$score="Polymorphism:$score";
					}elsif ($pred eq 'P'){
						$score="Polymorphism Automatic:$score";
					}elsif ($pred eq 'NA'){
						$score="NA:$score";
					}
				}elsif ($dbtype=~/ma/ && $keepline){
					if ($score<0.7925){
						$score="Neutral:$score";
					}elsif ($score<1.8925){
						$score="Low:$score";
					}elsif ($score<3.4925){
						$score="Medium:$score";
					}else{
						$score="High:$score";
					}
				}
			} elsif ($dbtype =~ m/^ljb23_/) {
				@record = split (/\t/, $_);my ($pred,$transformedValue);
				@record >= 5 or next;# die "Error: invalid record found in DB file $dbfile (at least 5 tab-delimited fields expected): <$_>\n";
				($chr, $start, $end, $ref, $obs, $score, $transformedValue, $pred) = @record;$pred||='';
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				if (defined $score and defined $score_threshold and not $keepline) {
					if ($reverse) {
						$score > $score_threshold and next;
					} else {
						$score < $score_threshold and next;
					}
				}
				if ($dbtype=~/pp2/ && $keepline){
					$pred=$transformedValue;
					if ($pred eq 'D'){
						$score="DAMAGING:$score";
					}elsif ($pred eq 'P'){
						$score="Probably DAMAGING:$score";
					}elsif ($pred eq 'B'){
						$score="Benign:$score";
					}elsif ($pred eq 'NA'){
						$score="NA:$score";
					}
				}elsif ($dbtype=~/mt/ && $keepline){
					# MutationTaster prediction, "A" ("disease_causing_automatic"), "D" ("disease_causing"), "N" ("polymorphism") or "P" ("polymorphism_automatic")
					if ($pred eq 'A'){
						$score="Disease Causing Automatic:$score";
					}elsif ($pred eq 'D'){
						$score="Disease Causing:$score";
					}elsif ($pred eq 'N'){
						$score="Polymorphism:$score";
					}elsif ($pred eq 'P'){
						$score="Polymorphism Automatic:$score";
					}elsif ($pred eq 'NA'){
						$score="NA:$score";
					}
				}elsif ($dbtype=~/ma/ && $keepline){
					if ($pred eq 'L'){
						$score="Low:$score";
					}elsif ($pred eq 'M'){
						$score="Medium:$score";
					}elsif ($pred eq 'N'){
						$score="Neutral:$score";
					}elsif ($pred eq 'H'){
						$score="High:$score";
					}elsif ($pred eq 'NA'){
						$score="NA:$score";
					}
				}
			} elsif ($dbtype =~ m/^ljb26_/) {
				@record = split (/\t/, $_);my ($pred,$transformedValue);
				@record >= 5 or next;# die "Error: invalid record found in DB file $dbfile (at least 5 tab-delimited fields expected): <$_>\n";
				($chr, $start, $end, $ref, $obs, $score, $pred) = @record;$pred||='';
				$chr=~s/chr//g;
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				if (defined $score and defined $score_threshold and not $keepline) {
					if ($reverse) {
						$score > $score_threshold and next;
					} else {
						$score < $score_threshold and next;
					}
				}
				# print "dbtype:$dbtype\n";
				if ($dbtype=~/(pp2|sift)/ && $keepline){
					if ($pred eq 'D'){
						$score="DAMAGING:$score";
					}elsif ($pred eq 'P'){
						$score="Probably DAMAGING:$score";
					}elsif ($pred eq 'B'){
						$score="Benign:$score";
					}elsif ($pred eq 'NA'){
						$score="NA:$score";
					}elsif($pred eq 'T'){
						$score="Tolerated:$score";
					}
				}elsif ($dbtype=~/mt/ && $keepline){
					# MutationTaster prediction, "A" ("disease_causing_automatic"), "D" ("disease_causing"), "N" ("polymorphism") or "P" ("polymorphism_automatic")
					if ($pred eq 'A'){
						$score="Disease Causing Automatic:$score";
					}elsif ($pred eq 'D'){
						$score="Disease Causing:$score";
					}elsif ($pred eq 'N'){
						$score="Polymorphism:$score";
					}elsif ($pred eq 'P'){
						$score="Polymorphism Automatic:$score";
					}elsif ($pred eq 'NA'){
						$score="NA:$score";
					}
				}elsif ($dbtype=~/ma/ && $keepline){
					if ($pred eq 'L'){
						$score="Low:$score";
					}elsif ($pred eq 'M'){
						$score="Medium:$score";
					}elsif ($pred eq 'N'){
						$score="Neutral:$score";
					}elsif ($pred eq 'H'){
						$score="High:$score";
					}elsif ($pred eq 'NA'){
						$score="NA:$score";
					}
				}elsif ($dbtype=~/(fathmm|MetaSVM|MetaLR)/i && $keepline){
					if ($pred eq 'T'){
						$score="Tolerated:$score";
					}elsif ($pred eq 'D'){
						$score="Deleterious:$score";
					}
				}elsif ($dbtype=~/cadd/ && $keepline){
					$score="$pred";	
				}
			}elsif ($dbtype=~/sift/){
				@record = split (/\t/, $_);
				@record >= 5 or next;#die "Error: invalid record found in DB file $dbfile (8 tab-delimited fields expected): <$_>\n";
				($chr, $start, $end, $ref, $obs, $score) = @record;
				if (defined $score and defined $score_threshold and not $keepline) {
					my ($prediction,$num,$median)=split(/[:\(]/,$score);
					if ($reverse) {
						$num < $score_threshold and next;
					} else {
						$num > $score_threshold and next;
					}
				}
				
				if ($chromosome) {
					$valichr{$chr} or next;
				}
			} elsif ($dbtype =~ m/^(rbpvar)/) {
				# die "Here";
				@record = split (/\t/, $_, -1);		#-1 is required before some dbSNP records have many empty tab fields in the end
				# print join "\n",@record;exit;
				$record[1] =~ s/^chr// or next;# "Error: invalid record found in DB file (2nd field should start with 'chr'): <$_>\n";
				($chr, $start, $end, $rsid, $strand, $ucscallele, $twoallele,$score) = @record[1,2,3,4,6,8,9,26];
				($record[13],$record[14])=(sprintf("%.02f",$record[13]),sprintf("%.02f",$record[14]));
				$start++ unless ($start==$end);#hmv added unless			#UCSC use zero-start system
				if ($chromosome) {
					$valichr{$chr} or  next;
				}
				my @allele = split (/\//, $twoallele);
				# @allele >= 2 or next;		#Jan 2011 modification
				if ($strand eq '-') {					#handle reverse strand annotation (the vast majority of records in dbSNP should be already in + strand)
					for my $i (0 .. @allele-1) {
						$allele[$i] = revcom ($allele[$i]);
					}
				}
				for my $i (0 .. @allele-1) {
					if ($ucscallele eq $allele[$i]) {
						@obs2 = @allele;
						splice (@obs2, $i, 1);
						for my $j (0 .. @obs2-1) {
							if ($keepline){
								push (@score2,"$score");
							}else{
								push @score2, $score;
							}
						}
					}else{
						next;
					}
				}
				# print "OBS?". $score2[0] . "\n";
				if (@obs2) {
					$obs = shift @obs2;
					$score = shift @score2;
				} else {
					$verbose and printerr ("Database error: wildtype base $ucscallele is not part of the allele description in <$_>\n");
					next;
				}
				# print "$score...\n";<STDIN>;
			} elsif ($dbtype =~ m/^(snp\d+|regulome\d+)/) {
				@record = split (/\t/, $_, -1);		#-1 is required before some dbSNP records have many empty tab fields in the end
				unshift(@record,'000') if ($dbtype=~/regulome141/);
				# print join "\n",@record;exit;
				if ($dbtype=~/snp/){
					@record == 18 or @record == 26 or @record==24 or  next;# "Error: invalid record found in dbSNP database file $dbfile (18 or 26 fields expected but found ${\(scalar @record)}): <$_>\n" . join("\n",@record);
				}
				$record[1] =~ s/^chr// or next;# "Error: invalid record found in DB file (2nd field should start with 'chr'): <$_>\n";
				($chr, $start, $end, $rsid, $strand, $ucscallele, $twoallele, $class) = @record[1,2,3,4,6,8,9,11];
				($record[13],$record[14])=(sprintf("%.02f",$record[13]),sprintf("%.02f",$record[14]));
				$start++ unless ($start==$end);#hmv added unless			#UCSC use zero-start system
				if ($chromosome) {

					$valichr{$chr} or  next;
				}
				unless ($class eq 'single' or $class eq 'deletion' or $class eq 'in-del' or $class eq 'insertion') {	#enum('unknown','single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion')
					# print "skipping $class($_)\n" and next;##changed this 8/17/15hv
					# next;
				}
		
				my @allele = split (/\//, $twoallele);
				
				#before Jan 2011, only di-allelic SNPs are handled in ANNOVAR
				#@allele == 2 or next;		#many entries have no allele information (for example, rs71010435)
				#in Jan 2011 version, I decided to handle tri-allelic and quad-allelic SNP as well
				
				@allele >= 2 or next;		#Jan 2011 modification
				if ($strand eq '-') {					#handle reverse strand annotation (the vast majority of records in dbSNP should be already in + strand)
					for my $i (0 .. @allele-1) {
						$allele[$i] = revcom ($allele[$i]);
					}
					#$ucscallele = revcom ($ucscallele);		#added Jan 24, 2011 (per Kevin Ha) removed Feb 10, 2011 (per Eric Stawiski)
					#note that some SNPs (e.g., rs28434453) may have multiple location in diferent chromosome or strand; I may want to handle this by a special flag in the future
					#585        chr1    13301   13302   rs28434453      0       -       C       C       C/T     genomic single etc...
					#1367    chr15   102517867       102517868       rs28434453      0       +       G       G       C/T     genomic single etc...
				}
				
				#in-del is usually annotated below, so they require special treatment
				#587     chr1    384538  384539  rs3971283       0       +       T       T       -/ATT   genomic in-del  unknown 0       0       unknown exact   3
				if ($class eq 'in-del') {					#indel are usually annotated as -/xxx, where xxx is the alternative allele
					$obs = length ($ucscallele) . $allele[1];		#prefix a number before the alleles, indicating block substitution
					defined $allele[1] or die "no allele 1 <$_>";
				} elsif ($class eq 'insertion') {
					$start--;
					$obs = "0$allele[1]";
				} elsif ($class eq 'deletion') {
					$obs = length ($ucscallele);
				} else {
					for my $i (0 .. @allele-1) {
						if ($ucscallele eq $allele[$i]) {
							@obs2 = @allele;
							splice (@obs2, $i, 1);
							for my $j (0 .. @obs2-1) {
								if ($keepline){
									if ($buildver!~/hg/){
										push (@score2,"$rsid");
									}else{
										push (@score2,"$rsid,$record[13]($record[14])");
									}
								}else{
									push @score2, $rsid;
								}
							}
						}
					}
					if (@obs2) {
						$obs = shift @obs2;
						$score = shift @score2;
					} else {
						$verbose and printerr ("Database error: wildtype base $ucscallele is not part of the allele description in <$_>\n");
						next;
					}
				}
				if ($maf_threshold and not $keepline) {
					if ($reverse) {					#the frequency is the non-reference allele frequency, which could exceed 0.5
						$record[13] > $maf_threshold and next;
					} else {
						$record[13] < $maf_threshold and next;
					}
				}
				$score = "$rsid";
				if ($keepline && $dbtype=~/snp/ && $buildver=~/hg/){
					$score.=":$record[13]($record[14])";
				}elsif($dbtype=~/regulome/){
					$score=$record[$#record];
				}
				# print "$score...\n";next;
			} elsif ($dbtype =~ m/(^1000g(\w+)|^nci60)/ or $dbtype =~ m/^1000g2010_(\w+)/ or $dbtype =~ m/^1000g201\d\w\w\w_(\w+)/ or $dbtype=~m/(AFR|AMR|EAS|EUR|SAS|ALL).sites.*/) {	#dbtype1 should NOT be used here
				@record = split (/\t/, $_);
				@record == 5 or @record == 6 or next;#die "Error: invalid record found in 1000G database file $dbfile (5 or 6 fields expected): <$_>\n";
				($chr, $start, $ref, $obs, $af) = @record;			#there is no "END" in 1000G input file
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				if ($maf_threshold and not $keepline) {
					if ($reverse) {					#the frequency is the non-reference allele frequency, which could exceed 0.5
						$af > $maf_threshold and next;
					} else {
						$af < $maf_threshold and next;
					}
				}
				$score = $af;
			}elsif($dbtype=~/exac0\d+(.*)$/){
				my $other=$1;
				@record = split (/\t/, $_);
				my @cols;
				($chr, $start, $end, $ref, $obs,@cols) = @record;
				my @exac_headers=(
					"ExAC_$other"."ALL",
					"ExAC_$other". "AFR",
					"ExAC_$other".'AMR',
					"ExAC_$other".'EAS',
					"ExAC_$other". 'FIN',
					"ExAC_$other". 'NFE',
					"ExAC_$other" . 'OTH',
					"ExAC_$other". 'SAS');
				if($colsWanted){
					if ($colsWanted=~/all/){
						$score='';
						for (my $j=0;$j<=$#exac_headers;$j++){
							$score.="$exac_headers[$j]=$cols[$j];";
						}
						chop $score if ($score=~/;$/);
					}else{
						my @multi=split(",",$colsWanted);
						$score='';
						for (my $j=0;$j<=$#multi;$j++){
							next if ($multi[$j]>7);#There are only 8 columns (0-based)
							$score.="$exac_headers[$multi[$j]]=$cols[$multi[$j]];";
						}
						chop $score if ($score=~/;$/);
					}
				}else{
					$score="$exac_headers[0]=$cols[0]";
				}
				
				$chr=~s/chr//g;
			}elsif ($dbtype=~/(clinvar|cosmicdb|cosmic70|ExAC)/i && !$force){
				@record = split (/\t/, $_);
				@record == 5 or @record == 6 or next;#die "Error: invalid record found in 1000G database file $dbfile (5 or 6 fields expected): <$_>\n";
				($chr, $start, $end, $ref, $obs, $score) = @record;			#there is no "END" in 1000G input file
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				$chr=~s/chr//g;
			}elsif ($dbtype=~/(icgc)/i){
				@record = split (/\t/, $_);
				@record == 5 or @record == 6 or next;#die "Error: invalid record found in 1000G database file $dbfile (5 or 6 fields expected): <$_>\n";
				my $icgc_id;
				($chr, $start, $end, $ref, $obs, $icgc_id, $score) = @record;			#there is no "END" in 1000G input file
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				$chr=~s/chr//g;
			}elsif ($dbtype=~/(cadd)/i){
				@record = split (/\t/, $_);
				@record == 5 or @record == 6 or next;#die "Error: invalid record found in 1000G database file $dbfile (5 or 6 fields expected): <$_>\n";
				my $phredScore;
				($chr, $start, $end, $ref, $obs, $phredScore, $score) = @record;			#there is no "END" in 1000G input file
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				$chr=~s/chr//g;
			}elsif ($dbtype=~/(kaviar|esp6500|provean|hapmap|cg69|hrcr1)/i){#20160301 changed
				@record = split (/\t/, $_);
				# @record == 5 or @record == 6 or next;#die "Error: invalid record found in 1000G database file $dbfile (5 or 6 fields expected): <$_>\n";
				my (@hrscore);
				($chr, $start, $end, $ref, $obs, $score,@hrscore) = @record;			#there is no "END" in 1000G input file
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				$chr=~s/chr//g;
			} elsif ($dbtype eq 'generic'||  ($dbtype=~/viz\d+/ && $force)) {#removed $dbtype=~/(cg69|hapmap|provean|esp6500|kaviar)/i || 20160301
				($chr, $start, $end, $ref, $obs, $score) = split (/\t/, uc $_);		#make sure to use upper case, as query is always in upper case
				$chr=~s/chr//i;
				defined $obs or next;#"Error: the generic database file must contains at least five tab-delimited fields per line (but observed line: $_)\n";
				defined $score or $score = "NA";
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				defined $obs or die "Error: invalid record found in DB file $dbfile (at least 5 fields expected for 'generic' dbtype): <$_>\n";
				if ($start == $end and $ref eq '-') {	#insertion
					$obs = "0$obs";
				}
				if ($obs eq '-') {			#deletion
					$obs = $end-$start+1;
				} elsif ($start != $end) {		#block substitution
					$obs = ($end-$start+1) . $obs;
				}
				if (defined $score and $score=~/^[\d\.]*$/ and defined $score_threshold and not $keepline) {
					if ($reverse) {
						$score > $score_threshold and next;
					} else {
						$score < $score_threshold and next;
					}
				}
			} elsif ($dbtype eq 'vcf') {			#vcf file is adopted by 1000 Genomes Project; it can describe both SNPs and indels, and it may contain both summary level statistics and individual level genotype calls
				($chr, $start, $rsid, $ref, $obs, $qual, $fil, $info) = split (/\t/, $_);
				if ($chromosome) {
					$valichr{$chr} or next;
				}
	
				my ($ac, $an);
				if ($info =~ m/AF=([\d\.]+)/) {
					($score, @score2) = split (/,/, $1);	#this is ideal case scenario; in reality, 1000G does not supply AF for alternative alleles
					if ($obs =~ m/(\w+),(\w+)/) {		#1000G November; this format is not really valid because it does not handle tri-allelic SNP
						($obs, @obs2) = split (/,/, $obs);
						if (@obs2 and not @score2) {
							@score2 = map {$score} @obs2;
						}
					}
				} elsif ($info =~ m/AC=(\S+?);AN=(\d+)/) {
					my ($alleles, $count) = ($1, $2);
					if ($alleles =~ m/^(\d+),(.+)/) {
						$score = sprintf ("%.3f", $1/$count);
						@score2 = split (/,/, $2);
						@score2 = map {sprintf("%.3f", $_/$count)} @score2;
						($obs, @obs2) = split (/,/, $obs);				#the obs is composed of two alleles
					} else {
						$af = sprintf ("%.3f", $alleles/$count);
						$score = $af;
						#this is an invalid record in 1000GJuly: 1       2266231 rs11589451      C       T,A     .       PASS    AA=c;AC=20;AN=120;DP=237
						if ($obs =~ m/(\w+),/) {
							$count_invalid_dbline++;
							$verbose and printerr "WARNING: Invalid input line found in $dbfile (more than one alleles are observed, but only one is annotated with allelic counts): <$_>\n";
							$obs = $1;		#2011jun: instead of using "next", I decided to still process this type of variants
						}
					}
				} else {
					$score = 'NA';
					if ($obs =~ m/^(\w+),/) {
						($obs, @obs2) = split (/,/, $obs);
						@score2 = map {'NA'} @obs2;
					}
				}
				
				($start, $ref, $obs) = reformatStartRefObs ($start, $ref, $obs, $_);
				
				for my $i (0 .. @obs2-1) {		#if there are tri-allelic variants or quad-allelic variants or more alleles
					($start, $ref, $obs2[$i]) = reformatStartRefObs ($start, $ref, $obs2[$i], $_);
				}
			}elsif($force){
				@record = split (/\t/, $_);
				my @cols;
				($chr, $start, $end, $ref, $obs,@cols) = @record;
				if($colsWanted){
					$score = '';
					my @multi=split(",",$colsWanted);
					if ($headersWanted){
						my @headersWanted = split(",",$headersWanted);
						if ($#headersWanted == $#multi){
							for (my $j=0;$j<=$#multi;$j++){
								$score.="$headersWanted[$j]=$cols[$multi[$j]];";
							}
						}else{
							for (my $j=0;$j<=$#multi;$j++){
								$score.="$cols[$multi[$j]];";
							}
						}
					}elsif ($colsWanted=~/all/){
						for (my $j=0;$j<=$#cols;$j++){
							$score.="$cols[$j];";
						}
					}else{
						for (my $j=0;$j<=$#multi;$j++){
							$score.="$cols[$multi[$j]];";
						}
					}
				}elsif($headersWanted && $colsWanted=~/all/){
					$score='';
					my @headersWanted = split(",",$headersWanted);
					for (my $j=0;$j<=$#cols;$j++){
						$score.="$headersWanted[$j]=$cols[$j];";
					}
				}else{
					$score="$cols[0]";
				}
				
				$chr=~s/chr//g;
			} else {
				die "invalid dbtype: $dbtype\n";
			}
			# print Dumper ($variant);exit;
			if ( $variant->{$chr, $start, $obs}) {
				my ($ref, @info) = split (/\n/, $variant->{$chr, $start, $obs});	#most likely, only one piece of information
				for my $i (0 .. @info-1) {
					my @lineitems=split(/\t/,$info[$i]);
					if ($keepline){
						$score=~s/[\r\n]//g;
						print DROPPED "$lineitems[0]\t$score\n";#line number and score
					}else{
						$info[$i]=~s/(LINE\d+:)//g;
						#join ("\t", $dbtype, $score), "\t",
						print DROPPED  join("\t",@lineitems[1..$#lineitems]), "\n";
						
					}
				}
				delete $variant->{$chr, $start, $obs};
				#delete $variant->{$chr, $start, $obs};
			}
			if (@obs2) {
				for my $j (0 .. @obs2-1) {
					if ($variant->{$chr, $start, $obs2[$j]}) {
						my ($ref, @info) = split (/\n/, $variant->{$chr, $start, $obs2[$j]});	#most likely, only one piece of information
						for my $i (0 .. @info-1) {
							my @lineitems=split(/\t/,$info[$i]);
							if ($keepline){
								print DROPPED "$lineitems[0]\t$score2[$j]\n";
							}else{
								$info[$i]=~s/(LINE\d+:)//g;
								#join ("\t", $dbtype, "$score2[$j]"), "\t",;
								print DROPPED join ("\t",@lineitems[1..$#lineitems]), "\n";
							}
						}
						delete $variant->{$chr, $start, $obs2[$j]};
					}
				}
			}
			$chunk_here += $line_length;
			if ( $chunk_here > $chunk_max ) {
				last;
			}
		} #end while (<DB>)
	} #end $bb
	#$variant = $outer;

	
	for my $key (sort keys %$variant) {
		my ($chr, $start, $obs) = split ($;, $key);		#hash key separator
		my ($ref, @info) = split (/\n/, $variant->{$key});
		my $len;
		if ($obs =~ m/^(\d+)(.*)/) {
			($len, $obs) = ($1, $2);
			$obs ||= '-';			#deletion
			if ($len) {
				$end = $start+$len-1;
			} else {
				$end = $start;
			}
		} else {
			$end = $start;
		}
		my $addon;
		for my $i (0 .. @info-1) {
			my @lineitems=split(/\t/,$info[$i]);
			if ($keepline){
				print FIL "$lineitems[0]\t-\n";
			}else{
				$info[$i]=~s/^(LINE){0,1}(\d+:{0,1})\t//g;
				print FIL "$info[$i]\n";
			}
		}
	}
	printerr "Done\n";
	$count_invalid_dbline and printerr "WARNING: $count_invalid_dbline lines in dbfile $dbfile were ignored due to invalid formats\n";
}

sub reformatStartRefObs {
	my ($start, $ref, $obs, $line) = @_;
	if (length ($ref) == 1 and length ($obs) == 1) {#single base substitution
		1;					#the obs and obs2 is already handled
	} elsif ($obs =~ m/^\-((\w)(\w*))$/) {		#deletion (1000G March)
		$2 eq $ref or $ref eq 'N' or die "Error: mismatch of deleted allele and reference allele: <$_>\n";
		$obs = length ($1);
	} elsif ($obs =~ m/^\+(\w+)$/) {		#insertion (1000G March)
		$obs = "0$1";
	} elsif ($ref =~ m/^[ACGTN]+$/ and $obs =~ m/^[ACGTN]+$/) {
		if (length ($ref) > length ($obs)) {		#deletion or block substitution with shorter size
			my $head = substr ($ref, 0, length ($obs));
			if ($head eq $obs) {
				$start += length ($obs);
				$obs = length ($ref) - length ($obs);
			} else {
				$obs = length ($ref) . $obs;
			}
		} else {					#insertion or block substitution
			my $head = substr ($obs, 0, length ($ref));
			if ($head eq $ref) {
				$start += (length ($ref)-1);
				$obs = '0' . substr ($obs, length ($ref));
			} else {
				$obs = length ($ref) . $obs;
			}
		}
	} else {
		die "Error: invalid record found in VCF file: <$line>\n";
	}
	return ($start, $ref, $obs);
}
	
sub annotateQueryByRegion {
	open (QUERY, $queryfile) or die "Error: cannot read from --queryfile ($queryfile): $!\n";
	open (OUT, ">$outfile.${buildver}_$dbtype1") or die "Error: cannot write to output file $outfile.${buildver}_$dbtype1: $!\n";
	open (INVALID, ">$outfile.invalid_input") or die "Error: cannot write to output file $outfile.invalid_input: $!\n";

	my ($regiondb, $parent) = ({}, {});
	
	if ($dbtype eq 'gff3') {
		($regiondb, $parent) = readGFF3RegionAnnotation ();
	} elsif ($dbtype eq 'bed') {
		($regiondb) = readBedRegionAnnotation ();
	} else {
		($regiondb) = readUCSCRegionAnnotation ();
	}
	my ($chr, $start, $end, $ref, $obs);
	my ($invalid);
	my ($linecount, $invalidcount) = qw/0 0/;
	
	$time and printerr "NOTICE: Current time (before examining variants) is ", scalar (localtime), "\n";
	print OUT "#$dbtype\n" if ($keepline and !$noheader);#hmv added
	my $printout='';#keeps track of the stuff to print..saves until buffer is full;
	while (<QUERY>) {
		s/[\r\n]+$//;
		
		if (m/^#/ and $comment) {				#comment line start with #, do not include this is $linecount
			$printout.= "#comment\t#comment\t$_\n" if (!$noheader);#added if statement 20131122 but printout existed
			next;
		}elsif(m/^#/){
			next;
		}
		$linecount++;
		if ($linecount%100000==0){print OUT $printout if ($printout);$printout='';}##hmv added
		$invalid = 0;						#reset invalid status

		my @nextline = split (/\s+/, $_);
		($chr, $start, $end, $ref, $obs) = @nextline[@avcolumn];
		if ( not (defined $chr and defined $start and defined $end and defined $ref and defined $obs)) {
			$invalid++;
		} else {
			($ref, $obs) = (uc $ref, uc $obs);
			$zerostart and $start++;
			$chr =~ s/^chr//i;
			if ($chr =~ m/[^\w]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {
				$invalid++;
			} elsif ($ref eq '-' and $obs eq '-' 		#both are empty allele
				or $ref =~ m/[^ACTG0\-]/ 		#non-standard nucleotide code
				or $obs =~ m/[^ACGT0\-]/ 		#non-standard nucleotide code
				or $start =~ m/[^\d]/ 			#start is not a number
				or $end =~ m/[^\d]/ 			#end is not a number
				or $start > $end			#start is more than end
				# or $ref ne '0' and $end-$start+1 != length ($ref) 	#length mismatch with ref
				or $ref eq '-' and $start != $end	#length mismatch for insertion
				) {		
				$invalid++;
			}
		}


		if ($invalid) {
			print INVALID $_, "\n";	#invalid record found
			$printout.="ERR\n" if ($keepline);
			$invalidcount++;
			next;
		}

		my $bin1 = int ($start/$genomebinsize);		#start bin
		my $bin2 = int ($end/$genomebinsize);		#end bin (usually same as start bin, unless the query is really big that spans multiple megabases)
		my ($foundhit, $score, $name);
		for my $bin ($bin1 .. $bin2) {
			for my $nextgene (@{$regiondb->{$chr, $bin}}) {
				my ($txstart, $txend, $txscore, $txname) = @$nextgene;
				#print "working on $txstart,$txend, $txscore,$txname,$chr for $start,$end...\n";<STDIN>;
				if ($end < $txstart) {
					#db:            <------------------------->
					#query: <--->
					last;						#if genomic region is too far away from end, end the search of the bins
				} elsif ($end <= $txend) {				#query contained completely within db region
					if ($start >= $txstart) {
						#db:      <-------------------------->
						#query:       <------------------>
					} else {					#query overlap but upstream of db region
						#db:       <------------------------->
						#query: <---------------------->
						if ($minqueryfrac) {
							if (($end-$txstart+1)/($end-$start+1) < $minqueryfrac) {
								next;
							}
						}
					}
					$foundhit++;
					$score ||= $txscore; $name ||= $txname;
					if ($score < $txscore) {
						$score = $txscore;
						$name=$txname;
					}
					if ($score == $txscore and defined $name and $name ne $txname) {
						$name .= ",$txname";
					}
					if ($dbtype1 eq 'cytoBand') {			#a new chromosome band is encountered
						$name ne $txname and $name .= ",$txname";
					}
				} elsif ($start <= $txend) {
					if ($start >= $txstart) {			#query overlap but downstream of db region
						#db:      <------------------------>
						#query:        <----------------------->
						if ($minqueryfrac) {
							if (($txend-$start+1)/($end-$start+1) < $minqueryfrac) {
								next;
							}
						}
					} else {
						#db region completely contained within query
						#db:      <------------------------->
						#query: <------------------------------>
						if ($minqueryfrac) {
							if (($txend-$txstart+1)/($end-$start+1) < $minqueryfrac) {
								next;
							}
						}
					}
					$foundhit++;
					$score ||= $txscore; $name ||= $txname;
					if ($score < $txscore) {
						$score = $txscore;
						$name=$txname;
					}
					if ($score == $txscore and defined $name and $name ne $txname) {
						$name .= ",$txname";
					}
					if ($dbtype1 eq 'cytoBand') {			#a new chromosome band is encountered
						$name ne $txname and $name .= ",$txname";
					}
				} else {
					#query            ---
					#gene  <-*----*->
				}
			}
		}
		$linecount =~ m/000000$/ and printerr "NOTICE: Finished processing $linecount variants in queryfile\n";
		#HMV added block
		my $prefix= my $ender='';#hv added, we do this so we concatenate results at the end if $keepline is set
		if (!$keepline){#this is for the pipeline view so we can concatenate across all databases
			$prefix="$dbtype\t";
			$ender="\t$_";
		}
		#####hmv added references to line and ender below
		if ($foundhit) {
			$name ||= '';
			my @name = split (/,/, $name);
			my %name = map {$_, 1} @name;
			@name = keys %name; 
			
			if ($dbtype1 eq 'cytoBand') {
				map {s/^chr//} @name;
				if (@name >= 2) {
					$name[$#name] =~ s/^\d+//;
					$name = $name[0] . '-' . $name[$#name];
				} else {
					$name = $name[0];
				}
				$printout.= "$prefix$name". "$ender\n";
			} else {
				$name = join (";", @name);
				#changed this 20130725 to remove "Name=" in the output
				$printout.= "$prefix". $score?"Score=$score;":"";$printout.= $name?"$name":"" ;$printout.= "$ender\n";
			}
		}elsif ($keepline){##hmv added elsif
			$printout.= "-\n";#print "blank\n";
		}else{
			print "$foundhit and ($keepline)\n";
		}
	}
	print OUT $printout if ($printout);
	close (QUERY);
	close (OUT);
	close (INVALID);
	
	$time and printerr "NOTICE: Current time (after examining variants) is ", scalar (localtime), "\n";
	
	printerr "NOTICE: Finished region-based annotation on $linecount genetic variants in $queryfile";
	if ($invalidcount) {
		printerr " (including $invalidcount with invalid format written to $outfile.invalid_input)";
	} else {
		unlink ("$outfile.invalid_input");
	}
	printerr "\n";
	printerr "NOTICE: Output files were written to $outfile.${buildver}_$dbtype1\n";
	printerr "ANNOVAR Completed:\n\t".scalar (localtime)."\n";
}
sub annotateQueryByRegionWithIdx {
	my $flag;
	my $dbfile = File::Spec->catfile ($dbloc, "${buildver}_$dbtype1.txt"); 
	-e "$dbfile.idx" or $flag=1;# "Error: required database $dbfile does not exists. Please use 'annotate_variation.pl -downdb $dbtype $dbloc' to download annotation database. LN".__LINE__."\n";
	if ($flag){
		system ("touch $dbfile.notIdx.ERR\n");
		annotateQueryByRegion();
		return;
	}
	open (QUERY, $queryfile) or die "Error: cannot read from --queryfile ($queryfile): $!\n";
	open (OUT, ">$outfile.${buildver}_$dbtype1") or die "Error: cannot write to output file $outfile.${buildver}_$dbtype1: $!\n";
	open (INVALID, ">$outfile.invalid_input") or die "Error: cannot write to output file $outfile.invalid_input: $!\n";
	my (%variant, $filedone, $batchdone);
	my ($linecount, $batchlinecount, $invalid, $invalidcount) = (0, 0);
	my ($chr, $start, $end, $ref, $obs, $info);
	print OUT "#$dbtype\n" if ($keepline && !$noheader);#hmv added#added noheader 20131122
	while (1) {
		$_ = <QUERY>;
		if (not defined $_) {
			$filedone++;
		} else {
			s/[\r\n]//g;
			
			if (m/^#/ and $comment) {				#comment line start with #, do not include this is $linecount
				print  "#$_\n" if ($keepline && !$noheader);#added 20131122
				next;
			}elsif(m/^#/){
				next;
			}else{
				$linecount++;
				
				$batchlinecount++;
				if ($batchlinecount == $batchsize) {
					$batchdone++;
				}
				
				if ($memfree or $memtotal) {		#if these arguments are specified
					if ($linecount =~ m/00000$/) {						#about 40Mb memory per 10k lines for a typical input dataset
						my ($availmem, $allmem) = currentAvailMemory();
						$verbose and printerr "NOTICE: Current available system memory is $availmem kb (this program uses $allmem bytes memory), after reading $linecount query\n";
						if ($availmem and $availmem <= $memfree+50_000) {		#some subsequent steps may take ~50Mb memory, so here we try to allocate some more memory
							$batchdone++;
						}
						if ($memtotal and $allmem >= $memtotal-50_000) {	#when --memtotal is specified, ensure that program use less memory
							$batchdone++;
						}
					}
				}
		
				$invalid = 0;						#reset invalid status
				my @nextline = split (/\s+/, $_);my @all_else;
				($chr, $start, $end, $ref, $obs,@all_else) = @nextline[@avcolumn];
				if ( not (defined $chr and defined $start and defined $end and defined $ref and defined $obs) and !$novar) {
					$invalid++;
				} else {
					($ref, $obs) = (uc $ref, uc $obs);
					$zerostart and $start++;
					$chr =~ s/^chr//i;
					if ($novar){
					}elsif ($chr =~ m/[^\w]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {
						$invalid++;
					} elsif ($ref eq '-' and $obs eq '-' 		#both are empty allele
						or $ref =~ m/[^ACTG0\-]/ 		#non-standard nucleotide code
						or $obs =~ m/[^ACGT0\-]/ 		#non-standard nucleotide code
						or $start =~ m/[^\d]/ 			#start is not a number
						or $end =~ m/[^\d]/ 			#end is not a number
						or $start > $end			#start is more than end
						or ($ref ne '0' and $end-$start+1 != length ($ref) )	#length mismatch with ref
						or ($ref eq '-' and $start != $end	)#length mismatch for insertion
						) {		
						$invalid++;
					}
				}				
				if ($invalid) {
					print INVALID $_, "\n";	#invalid record found
					$invalidcount++;
#					print OUT "ERR$invalid at".__LINE__."\n" if ($keepline);
					next;
				}
				if ($start == $end and $ref eq '-') {	#insertion
					$obs = "0$obs";
				} elsif ($obs eq '-') {			#deletion
					$obs = $end-$start+1;
				} elsif ($end>$start or $start==$end and length($obs)>1) {	#block substitution	#fixed the bug here 2011feb19
					$obs = ($end-$start+1) . $obs;
				}
				$chr=~s/chr//g;
				if (exists $variant{$chr, $start}) {
					$variant{$chr, $start } .= "\n$linecount\t".join ("\t",$chr,$start,$end,$ref,$obs,@all_else);
				} else {
					$variant{$chr, $start} = "$linecount\t".join ("\t",$chr,$start,$end,$ref,$obs,@all_else);
				}
			}
			
		}
		if ($filedone or $batchdone) {
			printerr "NOTICE: Processing next batch with ${\(scalar keys %variant)} unique variants in $batchlinecount input lines\n";
			annotateNextRegion (\%variant,$dbfile);
			%variant = ();
			$batchlinecount = 0;				#reset the line count for this batch
			$batchdone = 0;
		}
		if ($filedone) {
			last;
		}
	}
	close (INVALID); close (OUT); close (QUERY);
	$time and printerr "NOTICE: Current time (after examining variants) is ", scalar (localtime), "\n";
	
	printerr "NOTICE: Finished region-based annotation on $linecount genetic variants in $queryfile";
	if ($invalidcount) {
		printerr " (including $invalidcount with invalid format written to $outfile.invalid_input)";
	} else {
		unlink ("$outfile.invalid_input");
	}
	printerr "\n";
	printerr "NOTICE: Output files were written to $outfile.${buildver}_$dbtype1\n";
	printerr "ANNOVAR Completed:\n\t" . `date`;
	close LOG;
}
sub annotateNextRegion{
	my $variant=shift;
	my $dbfile=shift;
	my ($BIN, $DBSIZE) = (0, 0);
	my %index = ();my %regiondb;
	my $bb = {};			#a subset of %index, which corresponds to the input variants
	my $flag_idx_search = 0;	#indicate if index-based search algorithm is used (faster speed for a small number of input variants
	if ( -f "$dbfile.idx" ) {
		open(IDX, "$dbfile.idx") or die "Error: cannot read from input database index $dbfile.idx: $!\n";
		my $line = <IDX>;
		if (not $line =~ m/BIN\t(\d+)\t(\d+)/) {
			printerr "WARNING: Malformed database index file $dbfile.idx.\n";
		} elsif ($2 != -s $dbfile) {		#file version is different, do not use index file in this case
			printerr "WARNING: Your index file $dbfile.idx is out of date and will not be used. ANNOVAR can still generate correct results without index file.\n";
		} else {
			($BIN, $DBSIZE) = ($1, $2);
			my $last_chunk=0;#we do this to keep track of the last results
			my $last_chrom=0;
			while ( $line = <IDX> ) {
				$line =~ s/[\r\n]+$//;
				my ( $chrom, $pos, $offset0, $offset1 ) = split (/\t/, $line);
				defined $offset1 or next;				#invalid input line in the index file
				$chrom=~s/chr//i;#"chr$chrom" if ($chrom!~/chr/);
#				if ($chrom eq "$last_chrom"){
#					if ($last_chunk!=($pos-$BIN)){
#						$index{"$chrom\t$pos"} = $index{"$chrom\t$last_chunk"};
#					}else{
#						$index{"$chrom\t$pos"} = [$offset0, $offset1];
#					}
#				}else{
					# my $quart=int(($offset1-$offset0)/4);
					$index{"$chrom\t$pos"} = [ $offset0,$offset1];
#				}
#				$last_chunk=$pos;
				
			}
#			$index{"X$last_chrom"}=$last_chunk;#for the final chromosome
		}
		close (IDX);#print Dumper (\%index);<STDIN>;
		if (not %index) {
			printerr "WARNING: Unable to load database index successfully from file $dbfile.idx.\n";
		} else {
			foreach my $k ( sort keys %$variant ) {
				#my ( $chrom, $pos ) = split /\t/, $outer->{$k};
				#$chrom =~ s/.+\n//;
				my ($chrom, $pos) = split ($;, $k);
				$chrom=~s/^chr//i;
				my $bin = $pos - ( $pos % $BIN );
#				print "$k and $bin\n";
				
				if (!defined $index{"$chrom\t$bin"} ){
					my $fail=0; #keep track of failures so I don't keep looping
					until (defined  $index{"$chrom\t$bin"} || $bin<0 || $fail >10 ){
						$fail++;
						$bin-=$BIN;
					}
				}
#				print "$chrom,$pos,$bin at line ".__LINE__."\n"; 
				push (@{$bb->{ "$chrom\t$bin" }}, $variant->{ $k });
#				print Dumper ($variant);
#				print "DUMP:$chrom\t$bin:". Dumper ($bb->{ "$chrom\t$bin" })."\n|||$chrom\t$bin|||$k";<STDIN>;
			}
			
			if (scalar (keys %$bb) / scalar (keys %index) <= $indexfilter_threshold) {#HV2013_05
				$flag_idx_search++;
				printerr "NOTICE: Database index loaded. Total number of bins is ".  scalar (keys %index) . " and " . " number of bins to be scanned is " . scalar (keys %$bb) . "\n";
#				print STDERR Dumper ($bb);
			}
		}
	}else{
		system ("touch $dbfile.notIdx.ERR\n");
		annotateQueryByRegion();
		return;
	}
#	print Dumper (\%index);<STDIN>;
#	if (not $flag_idx_search) {
#		$bb = {1, [0, -s "$dbfile"]};
#		%index = (1, [0, -s "$dbfile"]);
#	}
#	print Dumper ($bb);die;
	my ($chr, $start, $end, $ref, $obs);
	my $exactPos=0;#for not really region annotations but for specific positions and occasionally a range ## does not check the alleles, too many variables to consider
	$exactPos=1 if ($dbtype=~/^avsift|generic|1000g_(ceu|yri|jptchb)|1000g2010_(ceu|yri|jptchb)|1000g20\d\d[a-z]{3}_[a-z]+|regulome|snp\d+|vcf|(ljb_\w+)$/i);
	my ($invalid);
	my ($linecount, $invalidcount) = qw/0 0/;
#	open (DB, $dbfile) or die "Cannot open $dbfile\n";
	printerr "NOTICE: Current time (before examining variants) is ", scalar (localtime), "\n";
	my %printout;
	my $abcc_aref=_getABCC_DBInfo($dbloc,$dbtype1);#use the lkup file to determine where the chr positions, name and score are

	# print Dumper (\%index);<STDIN>;
	foreach my $batch ( sort keys %$bb ) {
		# print "testing $batch| ".join ("\n",@{$$bb{$batch}})."\n";<STDIN>;
		next if ($#{$$bb{$batch}} == -1);
		#my ( $chrom, $bin ) = split /\t/, $b;
		my ( $chunk_min, $chunk_max );
		if (!exists $index{$batch}){
		} else{
			( $chunk_min, $chunk_max )= @{ $index{$batch} };

		# printerr "Processing bin=$batch min=$chunk_min max=$chunk_max\n";
		# print STDERR "$chunk_min,$batch\t" .Dumper $bb->{$batch};<STDIN>;
		#close(DB);
			open(DB, $dbfile);
			if ($chunk_min<0){$chunk_min=0;}
			my $success=seek( DB, $chunk_min, 0 );
			warn "could not seek at ($chunk_min) position". __LINE__. "\n" if (!$success);
			#$variant = $bb->{ $b };
			my $chunk_here = $chunk_min;
			%regiondb=();
			my ($chr_idx,$start_idx,$end_idx,$score_idx,$cutoff_idx);
			if ($$abcc_aref[2]=~/,/){
				($chr_idx,$start_idx,$end_idx)=split(",",$$abcc_aref[2]);
				printerr "There is something wrong with your index file and/or database; Expected $chr_idx,$start_idx, $end_idx for your chromosomal position\n" if ($chr_idx eq '' || !$start_idx || !$end_idx);
			}else{
				die "haven't coded for yet...$$abcc_aref[2]\n".join ("|",@{$abcc_aref})."\n";
			}
			if ($$abcc_aref[3]=~/,/){
				($score_idx,$cutoff_idx)=split(",",$$abcc_aref[3]);
			}else{
				$score_idx=$$abcc_aref[3] ;
			}
			while (<DB>) {
				chomp;
				next if (!$_);
				my @data=split(/\t/,$_);chomp $data[$#data];
				my ($feat_chr,$feat_start,$feat_end,$feat_name,$feat_score);
				($feat_chr,$feat_start,$feat_end)=($data[$chr_idx],$data[$start_idx],$data[$end_idx]);
				# print "\nworking on $feat_chr,$feat_start,$feat_end($chr_idx,$start_idx,$end_idx)\n";;<STDIN>;
				$feat_chr=~s/^chr//i; next if ($feat_chr!~/^[\dXYMT]{1,2}$/);
				next if ($feat_start!~/^\d+$/ || $feat_end!~/^\d+$/);
				$feat_start++ if ($dbfile=~/(snp|gwas)/ && $feat_start<$feat_end);#UCSC numbering 
				
				if ($score_idx){
					$feat_score="Score=$data[$score_idx];";
					defined $score_threshold and $data[$score_idx] < $score_threshold and $feat_score='';			#if --score_threshold is set, the low scoring segment will be skipped
					defined $normscore_threshold and $data[$cutoff_idx] < $normscore_threshold and $feat_score='';		#if --normscore_threshold is set, the low scoring segment will be skipped
				}
				if ($#data<$$abcc_aref[4]){next;}
				$feat_name=($$abcc_aref[4] && $data[$$abcc_aref[4]] ne '')?$data[$$abcc_aref[4]]:"No feat Name";
#				for my $i ( $feat_start .. $feat_end ) {
#					if ($i==$feat_start ){
#						$regiondb{"$feat_chr,$i"}.=($feat_score)?"$feat_name,$feat_score;": "$feat_name;";
#					}else{
#						$regiondb{"$feat_chr,$i\_"}.=($feat_score)?"$feat_name,$feat_score;": "$feat_name;";
#					}
#					last if $exactPos;
#				}
				if ($feat_start!=$feat_end  ){
					$regiondb{"$feat_chr,$feat_start,$feat_end"}.=($feat_score)?"$feat_name,$feat_score;": "$feat_name;";
				}else{
					$regiondb{"$feat_chr,$feat_start,$feat_end"}.=($feat_score)?"$feat_name,$feat_score;": "$feat_name;";
				}
				last if (tell(DB)>$chunk_max);
			}
			close DB;
		}
		my $found;
		foreach my $myline (@{$$bb{$batch}}){
			my @duplicates=split("\n",$myline);
			foreach my $varline (@duplicates){
				$found=0;#reset for each variant
				my ($qline,$qchr,$qstart,$qend,@other)=split ("\t",$varline);
#				for my $m ($qstart .. $qend ) {
				foreach my $regionAnnot (keys %regiondb){
					my ($r,$s,$e)=split(",",$regionAnnot);
					if ($s!~/^\d+$/ || $e!~/^\d+$/ || $qstart!~/^\d+$/ || $qend!~/^\d+$/){print join (",",$s,$e,$qstart,$qend). "at LN $myline($varline)"; <STDIN>;}
#					if (exists $regiondb{"$qchr,"} || exists $regiondb{"$qchr,$m\_"} ){
					if ( ($qstart>=$s && $qstart<=$e )|| ($qend>=$s &&$qend<=$e ) || ($qend<=$s && $qend>$e)){
						$found++;$printout{$qline}.=$regiondb{$regionAnnot};
					}
					last if ($found);
				}
				if (!$found){
					$printout{$qline}="-" if ($keepline);
					if ($qstart==177930835){print "not found ($qstart)";}
					if ($qstart==184005719){print "$batch...not found ($qstart)";print Dumper (\%regiondb);<STDIN>;}
				}elsif ($qstart==177930835){print " found ($qstart)";
				}elsif ($qstart==184005719){print " found ($qstart)";}
			}
		}
	}
	my $lastline='';
	foreach my $line (sort {$a <=> $b} keys %printout){
		if ($lastline eq ''){
		}elsif ($keepline && $lastline+1 !=$line){
			for my $j( $lastline+1 .. $line-1){
				print OUT "ERR\n" ;
			}
		}
		chop $printout{$line} if ($printout{$line}=~/[;,:]$/);
		print OUT "$printout{$line}\n";
		$lastline=$line;
	}
	return;
}
sub readGFF3RegionAnnotation {
	my ($dbfile);
	my ($regioncount, $dbcount) = (0, 0);
	my (@record, %regiondb, %parent);
	
	$dbfile = File::Spec->catfile ($dbloc, $gff3dbfile);
	-f $dbfile or die "Error: required database $dbfile does not exists. Please use 'annotate_variation.pl -downdb $dbtype $dbloc -buildver $buildver' to download annotation database.\n";
	
	open (DB, $dbfile) or die "Error: cannot read from database file $dbfile: $!\n";
	printerr "NOTICE: Reading annotation database $dbfile ... ";
	$_ = <DB>;
	$_ =~ m/^##gff-version\s+3/ or die "Error: invalid header line found in the GFF3 database $dbfile (expect to see '##gff-version 3'): <$_>\n";
	while (<DB>) {
		m/^#/ and next;			#skip comments line
		m/^##FASTA/ and last;		#reached the FASTA sequence section of GFF3 file
		$dbcount++;
		s/[\r\n]+$//;							#deleting the newline characters
		@record = split (/\t/, $_);
		@record == 9 or die "Error: invalid records found in the GFF3 database $dbfile (9 fields expected): <$_>\n";
		my ($chr, $start, $end, $score, $attribute) = @record[0,3,4,5,8]; 
		$chr=~s/^chr//i;			#sometimes the chr prefix is present and should be removed (query usually does not contain this chr prefix)
		my $name;
		defined $score_threshold and $score < $score_threshold and next;			#if --score_threshold is set, the low scoring segment will be skipped

		my @feature = split (/;/, $attribute);
		for my $i (0 .. @feature-1) {
			if ($feature[$i]=~/Name=(\S+)/){  ##DO NOT CHANGE THIS AS I NEED THIS FOR miRNA_ABCC 2013-05-31
				$name=$1;
			}elsif ($feature[$i]=~/ID=(\S+)/){
				$name=$1;
			}			
			#$feature[$i] =~ m/ID=(\S+)/ and $name = $1;
		}
		defined $name or die "Error: invalid record in GFF3 database $dbfile (ID field not found): <$_>\n";
		for my $i (0 .. @feature-1) {
			if ($feature[$i] =~ m/Parent=(.+)/) {
				my @parent = split (/,/, $1);
				for my $j (0 .. @parent-1) {
					$parent{$name} .= $parent[$j];
				}
			}
		}
		
		my ($bin1, $bin2) = (int($start/$genomebinsize), int($end/$genomebinsize));
		for my $nextbin ($bin1 .. $bin2) {
			push @{$regiondb{$chr, $nextbin}}, [$start, $end, $score, $name];
		}
		$regioncount++;
		if ($verbose and $dbcount =~ m/000000$/) {
			my ($availmem, $allmem) = currentAvailMemory();
			printerr "NOTICE: Current system available memory is $availmem kb (this ANNOVAR program used $allmem kb)\n";
		}
	}
	close (DB);
	for my $key (keys %regiondb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$regiondb{$key}} = sort {$a->[0] <=> $b->[0]} @{$regiondb{$key}};
	}
	printerr "Done with $regioncount regions\n";
	return (\%regiondb, \%parent);
}

sub readBedRegionAnnotation {
	my ($dbfile);
	my ($regioncount, $dbcount) = (0, 0);
	my (@record, %regiondb);
	my ($chr, $start, $end);
	
	$dbfile = File::Spec->catfile ($dbloc, $bedfile);

	-f $dbfile or die "Error: required bedfile $dbfile does not exists.\n";
	
	open (DB, $dbfile) or die "Error: cannot read from database file $dbfile: $!\n";
	printerr "NOTICE: Reading annotation database $dbfile ... ";

	while (<DB>) {
		$dbcount++;
		s/[\r\n]+$//;							#deleting the newline characters
		next if ($_=~/#/);
		@record = split (/\t/, $_);
		
		($chr, $start, $end) = @record;


		$chr =~ s/^chr//i;
		$start++;										#due to the zero-opening coordinate system in UCSC
		
		my ($bin1, $bin2) = (int($start/$genomebinsize), int($end/$genomebinsize));
		for my $nextbin ($bin1 .. $bin2) {
			push @{$regiondb{$chr, $nextbin}}, [$start, $end, 0, 'NA'];
		}
		$regioncount++;
		if ($verbose and $dbcount =~ m/000000$/) {
			my ($availmem, $allmem) = currentAvailMemory();
			printerr "NOTICE: Current system available memory is $availmem kb (this ANNOVAR program used $allmem kb)\n";
		}
	}
	close (DB);
	
	for my $key (keys %regiondb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$regiondb{$key}} = sort {$a->[0] <=> $b->[0]} @{$regiondb{$key}};
	}
	printerr "Done with $regioncount regions\n";
	return (\%regiondb);
}

sub readUCSCRegionAnnotation {
	my ($dbfile);
	my ($regioncount, $dbcount) = (0, 0);
	my (@record, %regiondb);
	my ($chr, $start, $end, $score, $normscore, $name);
	my (@positionCols, @scoreCols, @colsToOutput);
	
	if ($dbtype1 =~ m/^mce(\d+way)$/) {
		$dbfile = File::Spec->catfile ($dbloc, "${buildver}_phastConsElements$1.txt");
	} else {
		$dbfile = File::Spec->catfile ($dbloc, "${buildver}_$dbtype1.txt"); 
	}
	-f $dbfile or die "Error: required database $dbfile does not exists. Please use 'annotate_variation.pl -downdb $dbtype $dbloc' to download annotation database. LN".__LINE__."\n";
	my ($expectedLength,$pos,$colsOut);#=_getABCC_DBInfo($dbtype1);
#	#################$$$
#	### The following SWITCH structure is modified Jan 2011 to faciliate future expansion
#	### $expectedLength is the number of cols expected in each line
#	### @postionCols => location of ($chr,$start,$end) columns
#	### @scoreCols => location of ($score, $normscore) columns leave empty is set not present (then set to zero below) ; WARNING must be empty or of length 2
#	### @colsToOutPut => location of ($name) columns to put into $name concatinated with ":" below
#	if ($dbtype1 =~ m/^phastConsElements\d+way/) {
#		$expectedLength=6;
#		@positionCols=(1,2,3);
#		@scoreCols=(4,5);		#normalized score
#		@colsToOutput=(4);		#lod=xxx is the Name output
#	} elsif ($dbtype1 eq 'evofold') {
#		$expectedLength=10;
#		@positionCols=(1,2,3);
#		@scoreCols=(5,5);
#		@colsToOutput=(4);
#	} elsif ($dbtype1 eq 'tfbsConsSites') {
#		$expectedLength=8;
#		@positionCols=(1,2,3);
#		@scoreCols=(7,5);
#		@colsToOutput=(4);
#	} elsif ($dbtype1 eq 'wgRna') {
#		$expectedLength=10;
#		@positionCols=(1,2,3);
#		@scoreCols=(5,5);
#		@colsToOutput=(4);
#	} elsif ($dbtype1 eq 'targetScanS') {
#		$expectedLength=7;
#		@positionCols=(1,2,3);
#		@scoreCols=(5,5);
#		@colsToOutput=(4);
#	} elsif ($dbtype1 eq 'genomicSuperDups') {
#		$expectedLength=30;
#		@positionCols=(1,2,3);
#		@scoreCols=(27,27);
#		@colsToOutput=(4);
#	} elsif ($dbtype1 eq 'omimGene') {
#		$expectedLength=5;
#		@positionCols=(1,2,3);
#		@scoreCols=();
#		@colsToOutput=(4);
#	} elsif ($dbtype1 eq 'gwasCatalog') {
#		$expectedLength=23;
#		@positionCols=(1,2,3);
#		@scoreCols=();
#		@colsToOutput=(10);
#	} elsif ($dbtype1 eq 'dgv') {
#		$expectedLength=16;
#		@positionCols=(1,2,3);
#		@scoreCols=();
#		@colsToOutput=(4); 
#	} elsif ($dbtype1 eq 'cytoBand') {		#special handling required
#		$expectedLength=5;
#		@positionCols=(0,1,2);
#		@scoreCols=();
#		@colsToOutput=(0,3);
#	} elsif ($dbtype1 =~ m/^chr\w+_chainSelf$/) {	#example: chr1_selfChain
#		$expectedLength=13;
#		@positionCols=(2,4,5);
#		@scoreCols=(12,12);
#		@colsToOutput=(11);
#	} elsif ($dbtype1 =~ m/^chr\w+_chain\w+$/) {	#example: chr1_chainPanTro2
#		$expectedLength=12;
#		@positionCols=(2,4,5);
#		@scoreCols=();
#		@colsToOutput=(11);
#	} elsif ($dbtype1 =~ m/^snp\d+/) {	
#		#$expectedLength=18;		#dbSNP132 in hg19 has now 26 fields!
#		$expectedLength='';
#		@positionCols=(1,2,3);
#		@scoreCols=();
#		@colsToOutput=(4);
#	} elsif ($dbtype1 eq 'gerp++elem') {
#		$expectedLength = 7;
#		@positionCols=(0,1,2);
#		@scoreCols=();
#		@colsToOutput=(3);
#	} elsif ($dbtype1=~/(ALL|YRI|CEU|JPTCHB).sites/i){
#		$expectedLength = 5;
#		@positionCols=(0,1,1);
#		@scoreCols=();
#		@colsToOutput=(4);
#	}elsif ($dbtype1=~/sift/i){
#		$expectedLength='';
#		@positionCols=(0,1,2);
#		@scoreCols=(5);
#		@colsToOutput=(5);
#	} else {
#		#other UCSC format if file is not defined above
#		$expectedLength='';
#		@positionCols=(1,2,3);
#		@scoreCols=();
#		@colsToOutput=(4);
#	} 
	my $abcc_aref=_getABCC_DBInfo($dbloc,$dbtype);
	die "your database file does not exist $dbfile\n" if (ref($abcc_aref)!~/ARRAY/i);
	$expectedLength=$$abcc_aref[1];
	@positionCols=split(",",$$abcc_aref[2]);
	$scorecolumn=$$abcc_aref[3];
	@colsToOutput=split(",",$$abcc_aref[4]);
	if ($scorecolumn) {
		@scoreCols = split(",",$scorecolumn);
	}
	open (DB, $dbfile) or die "Error: cannot read from database file $dbfile: $!\n";
	printerr "NOTICE: Reading annotation database $dbfile ... ";
	my $line;my $last=0;
	while ($line= <DB>){
		next if ($line=~/^#/);
		$last=1 if $line!~/^#/;
		last if ($last);
	}
	if ($expectedLength eq '') {		# if DB is unknown "generic format" use first line to get $expectedLength : file rewound afterwards
		@record = split (/\t/, $line);
		$expectedLength=@record;
		seek (DB, 0, 0);
	};
	
	########$$ Check to see if user has defined columns to output (intergers or all allowed)
	if (defined $colsWanted) { 
		if ($colsWanted[0] eq 'all') {
			@colsToOutput= 0 .. ($expectedLength-1);
		} elsif ($colsWanted[0] eq 'none') {
			@colsToOutput = ();
		} else{
			@colsToOutput = @colsWanted;
		}
	};

	########$$ check that the columns requested exist in the current DB
	for my $i (0 .. @colsToOutput-1) {
		if ($colsToOutput[$i] > $expectedLength) {
			die "Error: The DB file $dbfile has only $expectedLength columns but output column $colsToOutput[$i] is requested by --colsWanted!\n";
		}
	}
	while (<DB>) {
		$dbcount++;
		s/[\r\n]+$//;							#deleting the newline characters
		next if ($_=~/^#/ || $_ eq '');
		@record = split (/\t/, $_, -1);					#-1 is required so that trailing empty strings are not discarded (this is common in dbSNP files)
		
		@record == $expectedLength or warn "Error: invalid record in dbfile $dbfile ($expectedLength fields expected): <$_>\n";
		#HMV inserted the following block and put into a for loop so that it would loop through multiple genomic coordiantes
		my @positionArr;
		if ($#positionCols<2){
			$record[$positionCols[0]]=~s/_random//g;#mm9:X_random:562627-562648:-
			if ($record[$positionCols[0]]!~/,/){#hmv inserted this is for TargetScan  mm9:9:119034023-119034036:+
				my $extra;
				($extra,$chr,$start,$end)=split(/[:-]/,$record[$positionCols[0]] );#this works for miRNAorg_targets_conserved_aug2010
				$end=~s/[\r\n]+$//g;
				if ($end=~/[\+]/ || $end eq ''){#this works for deepBase_miRDeep_miRNA 
					($chr,$start,$end)=split(/[:-]/,$record[$positionCols[0]] );#move everything down one
				}
				push(@positionArr,$chr,$start,$end);
			}elsif ($record[$positionCols[0]]=~/((chr){0,1}[\dXYUn]+):(\d+)\-(\d+),(\d+)\-(\d+)/i){#mm9:9:119034023-119034036,119036106-119036114:+
				$chr=$1;
				push(@positionArr,$chr,$3,$4,$chr,$5,$6);
			}else{
				die __LINE__."|$record[$positionCols[0]]\n";
			}
		}else{
			@positionArr=@record[@positionCols];
		}
		for (my $h=0;$h<$#positionArr;$h+=3){
		#end of hmv inserted block...delete for end bracket...
	#		($chr, $start, $end) = @record[@positionCols];
			($chr, $start, $end) = (@positionArr[$h..($h+2)]);
			if (@colsToOutput) {		#I think there should always be a Name in the output column
				$name = join (':', @record[@colsToOutput]);
			}
			
			if(@scoreCols){
				($score, $normscore)=(@record[@scoreCols])
			} else{
				($score, $normscore) = qw/0 0/;
			}
			
			#########$$ Unusual exceptions for phastCons
			if ($dbtype1 =~ m/^phastConsElements\d+way/ && $buildver=~/hg/) {
				$score =~ s/^lod=// or die "Error: invalid lod score designation (no 'lod=' found) in dbfile $dbfile: $score<$_>\n";
			} ##lod= in the score for conservation tracks
			
			#########$$ Unusual exceptions for cytoBand
			if ($dbtype1 eq 'cytoBand' and not defined $colsWanted) {	#the name for chromosome band is concatenated as single word
				$name =~ s/://;
			}elsif ($dbtype1=~/(ALL|YRI|CEU|JPTCHB).sites/){
				$name= $1.":$name";
			}
			
			defined $score_threshold and $score < $score_threshold and next;			#if --score_threshold is set, the low scoring segment will be skipped
			defined $normscore_threshold and $normscore < $normscore_threshold and next;		#if --normscore_threshold is set, the low scoring segment will be skipped
			
			$chr =~ s/^chr//i;
#			print "test:$start,$end\n";<STDIN>;
			$start++ unless ($start==$end);									#due to the zero-opening coordinate system in UCSC
			
			my ($bin1, $bin2) = (int($start/$genomebinsize), int($end/$genomebinsize));
			for my $nextbin ($bin1 .. $bin2) {
				if ($rawscore) {								#print out rawscore, rather than normalized score (default)
					$normscore = $score;
				}
				if (defined $name) {
					push @{$regiondb{$chr, $nextbin}}, [$start, $end, $normscore, $name];
				} else {			#name is not requested in the output
					push @{$regiondb{$chr, $nextbin}}, [$start, $end, $normscore];
				}
			}
			$regioncount++;
			if ($verbose and $dbcount =~ m/000000$/) {
				my ($availmem, $allmem) = currentAvailMemory();
				printerr "NOTICE: Current system available memory is $availmem kb (this ANNOVAR program used $allmem kb)\n";
			}
		}###hmv added for block
	}
	close (DB);
	
	for my $key (keys %regiondb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$regiondb{$key}} = sort {$a->[0] <=> $b->[0]} @{$regiondb{$key}};
	}
	printerr "Done with $regioncount regions";
	if (defined $score_threshold or $normscore_threshold) {
		printerr " (that passed --score_threhsold or --normscore_threshold from a total of $dbcount regions)\n";
	} else {
		printerr "\n";
	}
	return (\%regiondb);
}


sub translateDNA {
	my ($seq) = @_;
	my ($nt3, $protein);
	$seq = uc $seq;
	#length ($seq) % 3 == 0 or printerr "WARNING: length of DNA sequence to be translated is not multiples of 3: <length=${\(length $seq)}>\n";
	while ($seq =~ m/(...)/g) {
		defined $codon1{$1} or printerr "WARNING: invalid triplets found in DNA sequence to be translated: <$1>\n";
		$protein .= $codon1{$1};
	}
	return $protein;
}

sub translateRNA {
	my ($seq) = @_;
	my ($nt3, $protein);
	$seq = uc $seq;
	#length ($seq) % 3 == 0 or printerr "WARNING: length of RNA sequence to be translated is not multiples of 3: <length=${\(length $seq)}>\n";
	while ($seq =~ m/(...)/g) {
		defined $codonr1{$1} or printerr "WARNING: invalid triplets found in RNA sequence to be translated: <$1>\n";
		$protein .= $codonr1{$1};
	}
	return $protein;
}

sub revcom {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/acgtACGT/tgcaTGCA/;
	return ($seq);
}

sub readSeqFromFASTADB {
	my ($refseqvar) = @_;
	my (%seqhash);
	my $seqdbfile;
	
	#the four statements below should be condensed in the future (they are identical)
	$seqdbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1" . "Mrna.fa");

	my ($seqid, $curseq) = ('', '');

	-f $seqdbfile or die "Error: FASTA sequence file $seqdbfile does not exist. Please use 'annotate_variation.pl --downdb $dbtype $dbloc' download the database.\n";
	open (SEQ, $seqdbfile) or die "Error: cannot read from seqdbfile $seqdbfile: $!\n";
	printerr "NOTICE: Reading FASTA sequences from $seqdbfile ... ";
	while (<SEQ>) {
		if (m/^>(\S+)/) {
			if ($refseqvar->{$seqid}) {
				not defined $seqhash{$seqid} and $seqhash{$seqid} = $curseq;		#finish reading the sequence for seqid and save it (unless the sequence is already read from the file)
			}
			$seqid = $1;
			$curseq = '';
		} else {
			if ($refseqvar->{$seqid}) {
				s/[\r\n]+$//;
				$curseq .= uc $_;					#only use upper case characters
			}
		}
	}
	if ($refseqvar->{$seqid}) {			#finish the last sequence in the file
		not defined $seqhash{$seqid} and $seqhash{$seqid} = $curseq;
	}
	close (SEQ);
	printerr "Done with ", scalar keys %seqhash, " sequences\n";
	if (keys %seqhash < keys %$refseqvar) {
		my (@seqnotfound, @seqnotfound_example);
		for $seqid (keys %$refseqvar) {
			exists $seqhash{$seqid} or push @seqnotfound, $seqid;
		}
		printerr "WARNING: A total of ${\(scalar @seqnotfound)} sequences cannot be found in $seqdbfile";
		@seqnotfound_example = splice (@seqnotfound, 0, 3);
		printerr " (example: @seqnotfound_example)\n";
	}
	return (\%seqhash);
}

sub readKgXref {
	my ($inputfile) = @_;
	my (%gene_xref);
	open (XREF, $inputfile) or die "Error: cannot read from kgxref file $inputfile: $!\n";
	while (<XREF>) {
		m/^#/ and next;
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		@record == 8 or die "Error: invalid record found in knownGene cross-reference file (6 fields expected): <$_>\n";
		#some genes were given names that are prefixed with "Em:" which should be removed due to the presence of ":" in exonic variant annotation
		#Em:AC006547.7 Em:AC005003.4 Em:U62317.15 Em:AC008101.5 Em:AC004997.11 Em:U51561.2
		$record[4] =~ s/^Em:/Em./;
		if ($gene_xref{$record[0]}) {			#BC003168 occur twice in kgxref file (OSBPL10, BC003168)
			if ($gene_xref{$record[0]} =~ m/^(BC|AK)\d+$/) {
				$gene_xref{$record[0]} = $record[4];
			}
		} else {
			$gene_xref{$record[0]} = $record[4];
		}
	}
	close (XREF);
	return (\%gene_xref);
}

sub readUCSCGeneAnnotation {			#read RefGene annotation database from the UCSC Genome Browser, convert 0-based coordinates to 1-based coordinates
	my ($dbloc) = @_;
	my ($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes);
	my (%genedb, %geneidmap, %name2count, %cdslen, %mrnalen);
	my ($genecount, $ncgenecount) = (0, 0);
	
	my $dbfile;
	my $kgxref;
	my %iscoding;		#this gene name is a coding gene (if it has coding and noncoding transcripts, ignore all noncoding transcripts)
	
	if ($dbtype1 eq 'refGene') {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");
	} elsif ($dbtype1 eq 'knownGene') {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");
		my $kgxreffile = File::Spec->catfile($dbloc, $buildver . "_kgXref.txt");
		-f $kgxreffile or die "Error: the knownGene cross-reference file $kgxreffile does not exist. Please use 'annotate_variation.pl --downdb knownGene $dbloc' to download the database.\n";
		$kgxref = readKgXref ($kgxreffile);
	} elsif ($dbtype1 eq 'ensGene') {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");
	} else {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");		#added 2011feb18
		#die "FATAL ERROR: the dbype $dbtype1 is not supported in the readUCSCGeneAnnotation() subroutine.\n";		#commented 2011feb18
	}
	-f $dbfile or die "Error: The gene annotation database $dbfile does not exist. Please use 'annotate_variation.pl --downdb $dbtype $dbloc -build $buildver' to download the database.\n";

	open (GENEDB, $dbfile) or die "Error: cannot read from gene annotaion database $dbfile: $!\n";
	printerr "NOTICE: Reading gene annotation from $dbfile ... ";
	while (<GENEDB>) {
		s/[\r\n]+$//;							#deleting the newline characters
		my @record = split (/\t/, $_);

		if ($dbtype1 eq 'refGene') {
			@record == 16 or die "Error: invalid record in $dbfile (expecting 16 tab-delimited fields in refGene file): <$_>\n";
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes) = @record[1..15];		#human hg18, mouse
		} elsif ($dbtype1 eq 'knownGene') {
			@record >= 11 or die "Error: invalid record in $dbfile (>=11 fields expected in knownGene file): <$_>\n";	#mm8=11, hg18=hg19=12
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend) = @record[0..9];
			$name2 = $kgxref->{$name} || $name;
		} elsif ($dbtype1 eq 'ensGene') {
			@record == 16 or die "Error: invalid record in $dbfile (expecting 16 fields in ensGene file): <$_>\n";
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes) = @record[1..15];
		} else {
			@record >= 11 or die "Error: invalid record in $dbfile (>=11 fields expected in $dbtype1 gene definition file): <$_>\n";
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes) = @record[1..15];
			defined $name2 or $name2=$name;
			#die "FATAL ERROR: the --dbtype $dbtype is not supported in readUCSCGeneAnnotation() subroutine.\n";		#commented 2011feb18
		}
	
		#handle situations where the same transcript is mapped to several chromosomes or regions (for example, NM_019105 is mapped to chr6, chr6_cox_hap1, chr6_qbl_hap2; NM_002538 is mapped to chr5 positive and negative strand and also in chr5_h2_hap1)
		if ($chr =~ m/hap\d+$/) {
			next;			#this is a temporary solution on 2011feb19, to ignore alternative haplotype chromosomes
		}
	
		$chr =~ s/^chr//i or die "Error: invalid record found in $dbfile (chrom field not found): <$_>\n";						#UCSC always prefix "chr" to the chromosome identifier, so this is a good check to make sure that the file is the correct file
		$dbstrand eq '+' or $dbstrand eq '-' or die "Error: invalid dbstrand information found in $dbfile (dbstrand has to be + or -): <$_>\n";		#dbstrand is important to know and cannot be optional
		my @exonstart = split (/,/, $exonstart); 			#remove trailing comma
		my @exonend = split (/,/, $exonend);				#remove trailing comma
		$exoncount == @exonstart or die "Error: invalid record found in $dbfile (exoncount discordance): <$exoncount vs ${\(scalar @exonstart)}>\n";
		@exonstart == @exonend or die "Error: invalid record found in $dbfile (exonstart and exonend count discordance): <${\(scalar @exonstart)} vs ${\(scalar @exonend)}>\n";
		$txstart++; $cdsstart++; map {$_++} @exonstart;			#convert 0-based coordinate to 1-based coordinate

		#LOGIC here:
		#first calcluate mRNA length, and if the transcript maps to multiple locations with discordant mRNA length, only consider the leftmost chromosome and leftmost coordinate (because the FASTA file is sorted in this manner)

		my $cdslength = 0;
		my $mrnalength = 0;
		for my $i (0 .. @exonstart-1) {
			$mrnalength += $exonend[$i]-$exonstart[$i]+1;
		}
		for my $i (0 .. @exonstart-1) {					#this calculation is valid regardless of strand
			#$mrnalength += $exonend[$i]-$exonstart[$i]+1;
			if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
				if ($cdsend <= $exonend[$i]) {
					$cdslength = $cdsend-$cdsstart+1;
					last;
				} else {
					$cdslength += $exonend[$i]-$cdsstart+1;
					next;
				}
			}
			if ($cdslength and $cdsend < $exonstart[$i]) {
				die "FATAL ERROR: impossible scenario for $name in $dbfile (cdsend is less than exon start)";
			} elsif ($cdslength and $cdsend <= $exonend[$i]) {
				$cdslength += $cdsend-$exonstart[$i]+1;
				last;
			} elsif ($cdslength and $cdsend > $exonend[$i]) {
				$cdslength += $exonend[$i]-$exonstart[$i]+1;
			}
		}
		
		if ($cdsstart != $cdsend+1) {		#coding gene
			if (defined $mrnalen{$name} and $mrnalen{$name} != $mrnalength) {
				$verbose and printerr "WARNING: $name occurs more than once in $dbfile with different mRNA length. The first occurences with identical mRNA length will be uesd in analysis.\n";
				next;
			}
			
			
			if (defined $cdslen{$name} and $cdslen{$name} != $cdslength) {
				$verbose and printerr "WARNING: $name occurs more than once in $dbfile with different CDS length. The first occurences with identical CDS length will be uesd in analysis.\n";
				next;
			}
			
			$iscoding{$name2}++;		#name2 is a coding gene, and if there is a noncoding transcript, ignore such transcripts in future analysis
		} else {		#noncoding gene
			1;
		}
		
		$cdslen{$name} = $cdslength;
		$mrnalen{$name} = $mrnalength;
				
		my ($bin1, $bin2) = (int(($txstart - $neargene)/$genomebinsize), int(($txend + $neargene)/$genomebinsize));
		for my $nextbin ($bin1 .. $bin2) {
			push @{$genedb{$chr, $nextbin}}, [$name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, [@exonstart], [@exonend], $name2];
		}
		$geneidmap{$name} = $name2;
		$genecount++;
		$name2count{$name2}++;
		$cdsstart == $cdsend+1 and $ncgenecount++;			#non-coding gene has the same start and end site
	} 
	close (GENEDB);
	
	my %badgene;
	for my $key (keys %genedb) {
		my @newgenedb;
		for my $geneinfo (@{$genedb{$key}}) {
			if (not $cdslen{$geneinfo->[0]} and $iscoding{$geneinfo->[8]}) {
				$badgene{$geneinfo->[0]}++;
			} else {
				push @newgenedb, $geneinfo;
			}
		}
		@{$genedb{$key}} = @newgenedb;
	}
	
	for my $key (keys %genedb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$genedb{$key}} = sort {$a->[2] <=> $b->[2]} @{$genedb{$key}};
	}
	printerr "Done with $genecount transcripts (including $ncgenecount without coding sequence annotation) for ", scalar (keys %name2count), " unique genes\n";
	$verbose and printerr "NOTICE: ", scalar (keys %badgene), " noncoding transcripts will be ignored, because their associated genes have annotated coding transcript\n";
	return (\%genedb, \%geneidmap, \%cdslen, \%mrnalen);
}

sub downloadDB {
	my ($cwd, $msg, $sc);
	
	$cwd = Cwd::cwd();
	
	-w $dbloc or die "Error: the directory $dbloc is not writable by the current user\n";
	chdir ($dbloc) or die "Error: the directory $dbloc cannot be accessed\n";
	
	my (@urlin, @filein, @fileout, %fail);		#the fail hash contains index of files that fail to be downloaded
	my $count_success;
	my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
	if ($dbtype1 eq 'refGene') {
		push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/refGene.txt.gz";
		push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/refLink.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_refGeneMrna.fa.gz";
	} elsif ($dbtype1 eq 'knownGene') {
		push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/knownGene.txt.gz";
		push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/kgXref.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_knownGeneMrna.fa.gz";
	} elsif ($dbtype1  eq 'ensGene') {
		push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/ensGene.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_ensGeneMrna.fa.gz";
	} elsif ($dbtype1 eq 'seq') {
		push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/chromFa.zip";		#example: hg18, hg19
		push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/chromFa.tar.gz";	#example: panTro2
		push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/$buildver.fa.gz";	#example: bosTau4
	} elsif ($dbtype1 =~ m/^mce(\d+way)$/) {		#it could be 17 way, 28 way, 30 way, 44 way, etc, depending on genome and on build
		push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/phastConsElements$1.txt.gz";
	} elsif ($dbtype1 eq '1000g') {				#dbtype1 is same as queryfile, when --downdb is used
		$buildver eq 'hg18' or die "Error: currently the --dbtype of '1000g' only support --buildver of 'hg18'\n";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_CEU.sites.2009_04.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_YRI.sites.2009_04.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_JPTCHB.sites.2009_04.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_CEU.sites.2009_04.txt.idx.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_YRI.sites.2009_04.txt.idx.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_JPTCHB.sites.2009_04.txt.idx.gz";
		#push @urlin, "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/release/2009_04/CEU.sites.2009_04.gz";		#these lines are commented Sep 2011
		#push @urlin, "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/release/2009_04/YRI.sites.2009_04.gz";
		#push @urlin, "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/release/2009_04/JPTCHB.sites.2009_04.gz";
	} elsif ($dbtype1 eq '1000g2010') {			#dbtype1 is same as queryfile, when --downdb is used
		$buildver eq 'hg18' or die "Error: currently the --dbtype of '1000g2010' only support --buildver of 'hg18'\n";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_CEU.sites.2010_03.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_YRI.sites.2010_03.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_JPTCHB.sites.2010_03.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_CEU.sites.2010_03.txt.idx.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_YRI.sites.2010_03.txt.idx.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_JPTCHB.sites.2010_03.txt.idx.gz";
	} elsif ($dbtype1 eq '1000g2010jul') {			#dbtype1 is same as queryfile, when --downdb is used
		$buildver eq 'hg18' or die "Error: currently the --dbtype of '1000g2010jul' only support --buildver of 'hg18'\n";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_CEU.sites.2010_07.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_YRI.sites.2010_07.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_JPTCHB.sites.2010_07.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_CEU.sites.2010_07.txt.idx.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_YRI.sites.2010_07.txt.idx.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_JPTCHB.sites.2010_07.txt.idx.gz";
	} elsif ($dbtype1 eq '1000g2010nov') {
		$buildver eq 'hg19' or die "Error: currently the --dbtype of '1000g2010nov' only support --buildver of 'hg19'\n";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg19_ALL.sites.2010_11.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg19_ALL.sites.2010_11.txt.idx.gz";
	} elsif ($dbtype1 =~ m/^1000g(\d{4})(\w{3})$/) {
		$buildver eq 'hg19' or die "Error: currently the --dbtype of $dbtype1 only support --buildver of 'hg19'\n";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg19_ALL.sites.${1}_$monthhash{$2}.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg19_ALL.sites.${1}_$monthhash{$2}.txt.idx.gz";
	} elsif ($dbtype1 eq 'null') {
		1;
	} elsif ($dbtype1 eq 'avsift') {
		printerr "NOTICE: the --webfrom argument is set to 'annovar' automatically for -dbtype of 'avsift'\n";
		$buildver eq 'hg18' or $buildver eq 'hg19' or die "Error: currently the --dbtype of avsift only support --buildver of 'hg18' or 'hg19'\n";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_avsift.txt.gz";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_avsift.txt.idx.gz";
	} else {
		if ($webfrom) {
			if ($webfrom eq 'annovar') {
				push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_$dbtype1.txt.gz";
			} elsif ($webfrom eq 'ucsc') {
				push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/$dbtype1.txt.gz";
			} else {
				push @urlin, "$webfrom/$dbtype1.txt.gz";
			}
		} else {
			push @urlin, "ftp://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/$dbtype1.txt.gz";	#default goes to UCSC
		}
		if ($dbtype1 =~ m/^snp\d+$/ and $webfrom eq 'annovar' or $dbtype1 =~ m/^cg\d+$/ or $dbtype1 =~ m/^ljb_\w+$/) {			#load filter-index for dbSNP databases or CG databases
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_$dbtype1.txt.idx.gz";
		}
	}
	
	@filein = @urlin;
	map {s/.+\///} @filein;
	@fileout = @filein;
	map {s/\.gz$//; s/\.zip$//} @fileout;
	
	if ($wget) {
		$msg = qx/wget --help 2>&1/ || '';		#collect the output of the system command
	} else {
		$msg = '';					#when --nowget is specified, do not use wget to retrieve files from Internet
	}
	if ($msg =~ m/Usage/) {
		checkProgramUpdate ("wget");
		for my $i (0 .. @urlin-1) {
			printerr "NOTICE: Downloading annotation database $urlin[$i] ... ";
			if ($verbose) {
				$sc = "wget -t 1 -T 10 -O $filein[$i] $urlin[$i]";
			} else {
				$sc = "wget -t 1 -T 10 -q -O $filein[$i] $urlin[$i]";
			}
			if (system ($sc)) {	#time-out is 10 seconds, with 1 retry attempt
				printerr "Failed\n";
				$verbose and print "WARNING: unable to execute system command: <$sc>\n";
				unlink ($filein[$i]);		#delete the temporary files generated by wget
				$fail{$i}++;
			} else {
				printerr "OK\n";
				$count_success++;
			}
		}
	} else {
		eval {
			require Net::FTP;
			require LWP::UserAgent;
		};
		if ($@) {
			printerr "WARNING: cannot retrieve remote files automatically (by 'wget' command or by standard Net::FTP/LWP::UserAgent Perl module).\n";
			printerr "Please manually download the following file, uncompress the files to $dbloc directory, then add a ${buildver}_ prefix to the file names.\n";
			printerr join ("\n", @urlin), "\n";
			exit (100);
		}
		
		checkProgramUpdate ("lwp");
		my ($http, $ftp);
		for my $i (0 .. @urlin-1) {
			printerr "NOTICE: Downloading annotation database $urlin[$i] ... ";
			if ($urlin[$i] =~ m/^http/) {
				$http = LWP::UserAgent->new (timeout=>10, show_progress=>$verbose);
				$http->env_proxy;
				
				my $response = $http->get ($urlin[$i], ':content_file'=>$filein[$i]);
				if ($response->is_success) {
					printerr "Done\n";
					$count_success++;
				} else {
					printerr "Failed\n";
					$verbose and printerr "WARNING: cannot retrieve remote files ($urlin[$i]) via LWP::UserAgent Perl module: ", $response->status_line, "\n";
					$fail{$i}++;
				}
			} elsif ($urlin[$i] =~ m#^ftp://([^\\\/]+)#) {		#for hgdownload.cse.ucsc.edu, ftp-trace.ncbi.nih.gov, ftp.ensembl.org, etc
				my $urlroot = $1;
				if ($ftp = Net::FTP->new($urlroot, Timeout=>10, Debug=>$verbose)) {
					$ftp->login("anonymous", 'anonymous@');
					$ftp->binary();
					my $url = $urlin[$i];
					$url =~ s#ftp://[\w\.\-]+/##;		#remove the URL root
					if (not $ftp->get($url)) {
						printerr "Failed\n";
						$verbose and printerr "WARNING: cannot retrieve remote file ($url) in FTP server $urlroot\n";
						$fail{$i}++;
					} else {
						printerr "Done\n";
						$count_success++;
					}
				} else {
					printerr "Failed\n";
					$verbose and printerr "WARNING: cannot retrieve remote file ($urlin[$i]) via Net::FTP Perl module\n";
					$fail{$i}++;
				}
				
			} else {
				die "Error: The URL $urlin[$i] uses an unsupported protocol. Download cannot continue\n";
			}
		}
	}
	
	$count_success and printerr "NOTICE: Uncompressing downloaded files\n";
	for my $i (0 .. @filein-1) {
		$fail{$i} and next;
		if ($filein[$i] =~ m/\.zip$/) {
			$msg = qx/unzip --help 2>&1/ || '';		#collect the output of the system command
			if ($msg =~ m/Usage/i) {
				if ($verbose) {
					system ("unzip -o $filein[$i]");
				} else {
					system ("unzip -o -q $filein[$i]");
				}
			} else {
				printerr "ERROR: unzip is not installed in your system.\nPlease manually uncompress the files (@filein) at the $dbloc directory", $dbtype1 eq 'seq'?".\n":", and rename them by adding ${buildver}_ prefix to the file names.\n";
				exit (101);
			}
		} elsif ($filein[$i] =~ m/\.tar\.gz$/) {		#panTro2 FASTA sequence is stored as tar.gz rather than zip
			$msg = qx/tar --help 2>&1/ || '';		#collect the output of the system command
			if ($msg =~ m/Usage/i) {
				system ("tar -x -z -f $filein[$i]");
			} else {
				printerr "ERROR: tar/gunzip is not installed in your system.\nPlease manually uncompress the files (@filein) at the $dbloc directory", $dbtype1 eq 'seq'?".\n":", and rename them by adding ${buildver}_ prefix to the file names.\n";
				exit (102);
			}
		} elsif ($filein[$i] =~ m/\.gz$/) {
			$msg = qx/gunzip --help 2>&1/ || '';		#collect the output of the system command
			if ($msg =~ m/Usage/i) {
				system ("gunzip -f $filein[$i]");
			} else {
				printerr "ERROR: gunzip is not installed in your system.\nPlease manually uncompress the files (@filein) at the $dbloc directory", $dbtype1 eq 'seq'?".\n":", and rename them by adding ${buildver}_ prefix to the file names.\n";
				exit (103);
			}
		}
	}

	for my $i (0 .. @fileout-1) {
		$fail{$i} and next;				#skip the file that failed to be downloaded
		my $fileout = $fileout[$i];
		$dbtype1 eq 'seq' and next;			#the zip file contains dozens of FASTA files so cannot rename them automatically
		if (not $fileout =~ m/^${buildver}_/) {		#if the buildver is not the prefix of the files
			rename ($fileout, "${buildver}_$fileout") or die "Error: cannot rename $fileout to ${buildver}_$fileout\n";
			$fileout = "${buildver}_$fileout";
		}
		if (not $fileout =~ m/\.txt$/ and not $fileout =~ m/\.fa$/ and not $fileout =~ m/\.idx$/) {
			rename ($fileout, "$fileout.txt");
		}
	}
	
	$count_success and printerr "NOTICE: Finished downloading annotation files for $buildver build version, with files saved at the '$dbloc' directory\n";
	$cwd and chdir ($cwd);
	if (%fail) {
		my @failindex = keys %fail;
		if ($dbtype1 eq 'seq' and @failindex == 1) {	#not really a fail, because for seq, ANNOVAR attempts on tar.gz and zip file
			1;
		} else {
			printerr "WARNING: Some files cannot be downloaded, including ", join (', ', @urlin[@failindex]), "\n";
		}
		
		for my $index (@failindex) {
			if ($urlin[$index] =~ m#^http://www\.openbioinformatics\.org.+Mrna.fa.gz$#) {
				printerr "---------------------------ADDITIONAL PROCEDURE---------------------------\n";
				printerr "--------------------------------------------------------------------------\n";
				printerr "NOTICE: the FASTA file $urlin[$index] is not available to download but can be generated by the ANNOVAR software. ";
				printerr "PLEASE RUN THE FOLLOWING TWO COMMANDS CONSECUTIVELY TO GENERATE THE FASTA FILES:\n\n";
				printerr "\tannotate_variation.pl --buildver $buildver --downdb seq $dbloc/${buildver}_seq\n";
				printerr "\tretrieve_seq_from_fasta.pl $dbloc/${buildver}_$dbtype1.txt -seqdir $dbloc/${buildver}_seq -format $dbtype1 -outfile $dbloc/${buildver}_${dbtype1}Mrna.fa\n";
				printerr "--------------------------------------------------------------------------\n";
				printerr "--------------------------------------------------------------------------\n";
			}
		}
	}	
}

sub currentAvailMemory {
	my ($availmem, $allmem) = (0, 0);
	if ($^O eq "MSWin32") {		#no easy solution to get available memory from Windows.
		($availmem, $allmem) = (0, 0);
	} elsif ($^O eq 'linux' or $^O eq 'aix' or $^O eq 'solaris') {
		if (open (TOP, "top -b -n 1 2>&1 |")) {
			my $index;
			while (<TOP>) {
				if (m/^Mem:.+\s(\d+)k free/) {
					$availmem = $1;
				}
				s/^\s+//;
				my @field = split (/\s+/, $_);
				@field >= 10 or next;			#make sure that the PID lines are reached
				if ($field[0] eq 'PID') {
					for my $i (0 .. @field-1) {
						$field[$i] eq 'RES' and $index = $i;
					}
				}
				if ($field[0] eq $$) {
					defined $index or die "Error: invalid output from top command: the line with PID and RES is not found\n";
					$allmem = $field[$index];
					if ($allmem =~ m/^([\d\.]+)(\w)$/) {
						if ($2 eq 'g') {
							$allmem = $1 * 1_000_000;
						} elsif ($2 eq 'm') {
							$allmem = $1 * 1_000;
						} elsif ($2 eq 'k') {
							$allmem = $1;
						} else {
							printerr "WARNING: unrecognizable output from top command: <$_>\n";
						}
					}
					last;
				}
			}
		}
	} else {
		($availmem, $allmem) = (0, 0);
	}
	return ($availmem, $allmem);
}

sub printerr{
	if ($silent){
		print STDERR @_;
	}
	print LOG @_;
}

sub checkProgramUpdate {
	my ($method) = @_;
	my $sc;
	my ($curdate, $webdate, $webdate1) = $LAST_CHANGED_DATE;
	my (@webcontent);
	$method eq 'wget' or $method eq 'lwp' or die "Error: update checking method can be only 'wget' or 'lwp'";
	printerr "NOTICE: Web-based checking to see whether ANNOVAR new version is available ... ";
	$LAST_CHANGED_DATE =~ m/LastChangedDate: (\d+)\-(\d+)-(\d+)/ or printerr "Failed\n" and return; 
	$curdate = $1.$2.$3;
	if ($method eq 'wget') {
		$sc = "wget -t 1 -T 10 -q -O .annovar_date http://www.openbioinformatics.org/annovar/download/annovar_date";
		if (system ($sc)) {
			printerr "Failed\n";
			return;
		} else {
			if (not open (AVDATE, ".annovar_date")) {
				printerr "Cannot access version information\n";
			} else {
				printerr "Done\n";
				@webcontent = <AVDATE>;		#$LAST_CHANGED_DATE =	'$LastChangedDate: 2011-10-05 09:10:39 -0700 (Wed, 05 Oct 2011) $';
				close (AVDATE);
				unlink (".annovar_date");
			}
		}
	} elsif ($method eq 'lwp') {
		my $http = LWP::UserAgent->new (timeout=>10);
		$http->env_proxy;
		my $response = $http->get("http://www.openbioinformatics.org/annovar/download/annovar_date");
		if ($response->is_success) {
			printerr "Done\n";
			$_ = $response->decoded_content;
			@webcontent = split (/\n/, $_);
		} else {
			printerr "Failed\n";
			return;
		}
	}
	
	$webdate = $webcontent[0];
	$webdate =~ s/[\r\n]+$//;
	$webdate1 = $webdate;
	$webdate1 =~ s/\-//g;			#remove the - sign in webdate
	if ($curdate < $webdate1) {
		printerr "----------------------------UPDATE AVAILABLE------------------------------\n";
		printerr "--------------------------------------------------------------------------\n";
		printerr "WARNING: A new version of ANNOVAR (dated $webdate) is available!\n";
		printerr "         Download from http://www.openbioinformatics.org/annovar/\n";
	
		if (@webcontent >= 2) {
			printerr "Changes made in the $webdate version:\n";
			for my $i (1 .. @webcontent-1) {
				if ($webcontent[$i] =~ m/^(\d{4})\-(\d{2})\-(\d{2})[\r\n]+$/) {
					$webdate = "$1-$2-$3";
					$webdate1 = "$1$2$3";
					if ($curdate >= $webdate1) {	#the current version is more recent than this date
						last;
					} else {
						printerr "Changes made in the $webdate version:\n";
					}
				} else {
					printerr "         * $webcontent[$i]";
				}
			}
		}
		printerr "--------------------------------------------------------------------------\n";
		printerr "--------------------------------------------------------------------------\n";
	}
}

=head1 SYNOPSIS 

 annotate_variation_ABCC.pl [arguments] <query-file|table-name> <database-location>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
        
        Arguments to download databases or perform annotations
            --downdb			download UCSC Genome Browser annotation database
            --geneanno			annotate variants by functional consequences on genes
            --regionanno		annotate variants by targetting specific genomics regions
            --filter			filter variants based on a position list
            --webfrom <string>		specify the source of database (default usually works fine)
        
        Arguments to control input and output
            --outfile <file>		output file prefix
            --zerostart			input query file uses half-open zero-start coordinate
            --dbtype <string>		database type
            --buildver <string>		genome build version (default: hg18 for human)
            --gff3dbfile <file>		specify the GFF3 DB file used in region-based annotation
            --genericdbfile <file>	specify the generic DB file used in filter-based annotation
            --vcfdbfile <file>		specify the DB file in VCF format in filter-based annotation
            --bedfile <file>		specify a BED file in region-based annotation
            --time			print out local time during program run
            --separate			separately print out all function of a variant (default: one line per variant)
            --colsWanted <string>	specify which columns to output in -regionanno by comma-delimited numbers
            --comment			print out comment line (those starting with #) in output files 
            --scorecolumn <int>		the column with scores in database file (for region-based annotation)
            --exonsort			sort the exon number in output line (for gene-based annotation)
            --transcript_function	use transcript name rather than gene name in gene-based annotation output
            --hgvs			use HGVS format for exonic annotation (c.122C>T rather than c.C122T)
        
        Arguments to fine-tune the annotation procedure
            --batchsize <int>		batch size for processing variants per batch (default: 5m)
            --genomebinsize <int>	bin size to speed up search (default: 100k for -geneanno, 10k for -regionanno)
            --expandbin <int>		check nearby bin to find neighboring genes (default: 2m/genomebinsize)
            --neargene <int>		distance threshold to define upstream/downstream of a gene
            --score_threshold <float>	minimum score of DB regions to use in annotation
            --reverse			reverse directionality to compare to score_threshold
            --normscore_threshold <float> minimum normalized score of DB regions to use in annotation
            --rawscore			output includes the raw score (not normalized score) in UCSC Browser Track
            --minqueryfrac <float>	minimum percentage of query overlap to define match to DB (default: 0)
            --splicing_threshold <int>	distance between splicing variants and exon/intron boundary (default: 2)
            --maf_threshold <float>	filter 1000G variants with MAF above this threshold (default: 0)
            --sift_threshold <float>	SIFT threshold for deleterious prediction (default: 0.05)
            --precedence <string>	comma-delimited to specify precedence of variant function (default: exonic>intronic...)
       
       Arguments to control memory usage
            --memfree <int>		ensure minimum amount of free system memory (default: 100000, in the order of kb)
            --memtotal <int>		limit total amount of memory used by ANNOVAR (default: 0, unlimited, in the order of kb)
            --chromosome <string>	examine these specific chromosomes in database file
            

 Function: annotate a list of genetic variants against genome annotation 
 databases saved at local disk.
 
 Example: #download gene annotation database (for hg18 build) and save to humandb/ directory
 	  annotate_variation.pl -downdb gene humandb/
 	  annotate_variation.pl -buildver mm9 -downdb mce30way mousedb/
 	  annotate_variation.pl -downdb snp130 humandb/
 	
 	  #gene-based annotation of variants in the varlist file (by default --geneanno is ON)
 	  annotate_variation.pl ex1.human humandb/
 	  
 	  #region-based annotate variants
 	  annotate_variation.pl -regionanno -dbtype mce44way ex1.human humandb/
 	  annotate_variation.pl -regionanno -dbtype gff3 -gff3dbfile tfbs.gff3 ex1.human humandb/
 	  
 	  #filter rare or unreported variants (in 1000G/dbSNP) or predicted deleterious variants
 	  annotate_variation.pl -filter -dbtype 1000g_ceu -maf 0.01 ex1.human humandb/
 	  annotate_variation.pl -filter -dbtype snp130 ex1.human humandb/
 	  annotate_variation.pl -filter -dbtype avsift ex1.human humandb/
 
 Version: $LastChangedDate: 2011-10-05 09:10:39 -0700 (Wed, 05 Oct 2011) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--downdb>

download annotation databases from UCSC Genome Browser, Ensembl, 1000 Genomes 
Project or other resources. The annotation files in this database are required 
for the functional annotation of variants.

=item B<--geneanno>

perform gene-based annotation. For each variant, examine whether it hit exon, 
intron, intergenic region, or close to a transcript, or hit a non-coding RNA 
gene, or is located in a untranslated region. In addition, for an exonic variant, 
determine whether it causes splicing change, non-synonymous amino acid change, 
synonymous amino acid change or frameshift changes.

=item B<--regionanno>

perform region-based annotation. For each variant, examine whether it overlaps 
with a specific genomic region, such as the most conserved elements, the 
predicted transcription factor binding sites, the specific cytogeneic bands, the 
evolutionarily conserved RNA secondary structures and so on.

=item B<--filter>

perform filter-based annotation. For each variants, filter it against a 
variation database, such as the 1000 Genomes Project database and the dbSNP 
database, and identify a subset that have not been reported in these databases 
as novel variants.

=item B<--outfile>

specify the output file prefix. Several output files will be generated using 
this prefix and different suffixes. A directory name can also be specified as 
part of the argument, so that the output files can be written to a different 
directory than the current directory.

=item B<--zerostart>

utilize the half-open zero-start coordinate system that is used by many UCSC 
Genome Browser annotation tables. By default, the 1-based coordinate system will 
be used.

=item B<--dbtype>

specify the database type to be used in gene-based, region-based or filter-based 
annotations. For gene-based annotation, by default refGene annotations from the 
UCSC Genome Browser will be used for annotating variants. However, users can 
switch to utilize Ensembl annotations instead, or use the UCSC Gene annotations. 
In general, RefSeq gene annotations are more conservative, and UCSC Gene 
annotations are most liberal with many predicted genes and non-coding RNAs. For 
region-based annotations, users can select any UCSC annotation databases (by 
providing the database name), or alternatively select a Generic Feature Format 
version 3 (GFF3) formatted file for annotation (by providing 'gff3' as the --
dbtype and providing the --gff3dbfile argument). For filter-based annotations, 
users can select a dbSNP file, a 1000G file, a generic format file (with simple 
columns including chr, start, end, reference, observed, score), a VCF format 
(which is the current popular format for variants exchange), or a avsift format 
(which is identital to the generic format but is provided for convenience).

=item B<--buildver>

genome build version to use. By default, the hg18 build for human genome is 
used. The build version will be used by ANNOVAR to identify corresponding database files 
automatically, for example, when gene-based annotation is used for hg18 build, 
ANNOVAR will search for the hg18_refGene.txt file, but if the hg19 is used as --
buildver, ANNOVAR will examine hg19_refGene.txt instead.

=item B<--gff3dbfile>

specify the GFF3-formatted database file used in the region-based annotation.

=item B<--genericdbfile>

specify the generic format database file used in the filter-based annotation.

=item B<--vcfdbfile>

specify the database file in VCF format in the filter-based annotation. VCF has 
been a popular format for summarizing SNP and indel calls in a population of 
samples, and has been adopted by 1000 Genomes Project in their most recent data 
release.

=item B<--time>

print out the local time during execution of the program

=item B<--separate>

for gene-based annotation, separate the effects of each variant, so that each 
effect (intronic, exonic, splicing) is printed in one output line. By default, 
all effects are printed in the same line, in the comma-separated form of 
'UTR3,UTR5' or 'exonic,splicing'.

=item B<--colsWanted>

specify which columns are desired in the output for -regionanno. By default, 
ANNOVAR inteligently selects the columns based on the DB type. However, users 
can use a list of comma-delimited numbers, or use 'all', or use 'none', to 
request custom output columns.

=item B<--comment>

specify that the program should include comment lines in the output files. 
Comment lines are defined as any line starting with #. By default, these lines 
are not recognized as valid ANNOVAR input and are therefore written to the 
INVALID_INPUT file. This argument can be very useful to keep columns headers in 
the output file, if the input file use comment line to flag the column headers 
(usually the first line in the input file).

=item B<--scorecolumn>

specify the the column with desired output scores in UCSC database file (for 
region-based annotation). The default usually works okay.

=item B<--exonsort>

sort the exon number in output line in the exonic_variant_function file during 
gene-based annotation. If a mutation affects multiple transcripts, the ones with 
the smaller exon number will be printed before the transcript with larger exon 
number in the output.

=item B<--batchsize>

this argument specifies the batch size for processing variants by gene-based 
annotation. Normally 5 million variants (usually one human genome will have 
about 3-5 million variants depending on ethnicity) are annotated as a batch, to 
reduce the amounts of memory. The users can adjust the parameters: larger values 
make the program slightly faster, at the expense of slightly larger memory 
requirements. In a 64bit computer, the default settings usually take 1GB memory 
for gene-based annotation for human genome for a typical query file, but this 
depends on the complexity of the query (note that the query has a few required 
fields, but may have many optional fields and those fields need to be read and 
kept in memory).

=item B<--genomebinsize>

the bin size of genome to speed up search. By default 100kb is used for gene-
based annotation, so that variant annotation focused on specific bins only 
(based on the start-end site of a given variant), rather than searching the 
entire chromosomes for each variant. By default 10kb is used for region-based 
annotation. The filter-based annotations look for variants directly so no bin is 
used.

=item B<--expandbin>

expand bin to both sides to find neighboring genes/regions. For gene-based 
annotation, ANNOVAR tries to find nearby genes for any intergenic variant, with 
a maximum number of nearby bins to search. By default, ANNOVAR will 
automatically set this argument to search 2 megabases to the left and right of 
the variant in genome.

=item B<--neargene>

the distance threshold to define whether a variant is in the upstream or 
downstream region of a gene. By default 1 kilobase from the start or end site of 
a transcript is defined as upstream or downstream, respectively. This is useful, 
for example, when one wants to identify variants that are located in the 
promoter regions of genes across the genome.

=item B<--score_threshold>

the minimum score to consider when examining region-based annotations on UCSC 
Genome Browser tables. Some tables do not have such scores and this argument 
will not be effective.

=item B<--normscore_threshold>

the minimum normalized score to consider when examining region-based annotations 
on UCSC Genome Browser tables. The normalized score is calculated by UCSC, 
ranging from 0 to 1000, to make visualization easier. Some tables do not have 
such scores and this argument will not be effective.

=item B<--rawscore>

for region-based annotation, print out raw scores from UCSC Genome Browser 
tables, rather than normalized scores. By default, normalized scores are printed 
in the output files. Normalized scores are compiled by UCSC Genome Browser for 
each track, and they usually range from 0 to 1000, but there are some 
exceptions.

=item B<--minqueryfrac>

The minimum fraction of overlap between a query and a database record to decide 
on their match. By default, any overlap is regarded as a match, but this may not 
work best when query consist of large copy number variants.

=item B<--splicing_threshold>

distance between splicing variants and exon/intron boundary, to claim that a 
variant is a splicing variant. By default, 2bp is used. ANNOVAR is relatively 
more stringent than some other software to claim variant as regulating splicing. 
In addition, if a variant is an exonic variant, it will not be reported as 
splicing variant even if it is within 2bp to an exon/intron boundary.

=item B<--maf_threshold>

the minor allele frequency (MAF) threshold to be used in the filter-based 
annotation for the 1000 Genomes Project databases. By default, any variant 
annotated in the 1000G will be used in filtering.

=item B<--memfree>

the minimum amount of free system memory that ANNOVAR should ensure to have. By 
default, if ANNOVAR takes too much memory such that only 100Mb system memory is 
available, ANNOVAR will stop reading annotation database into memory, and will 
start annotation procedure, and then clear the memory, and then read the next 
block of annotation database again. This argument ensures that ANNOVAR will not 
attempt to use virtual memory in the system (which makes ANNOVAR extremely 
slow).

=item B<--memtotal>

the total amount of memory that ANNOVAR should use at most. By default, this 
value is zero, meaning that there is no limit on that. Decreasing this threshold 
reduce the memory requirement by ANNOVAR, but may increase the execution time.

=item B<--chromosome>

examine these specific chromosomes in database file. The argument takes comma-
delimited values, and the dash can be correctly recognized. For example, 5-10,X 
represent chromosome 5 through chromosome 10 plus chromosome X.

=back

=head1 DESCRIPTION

ANNOVAR is a software tool that can be used to functionally annotate a list of genetic variants, 
possibly generated from next-generation sequencing experiments. For example, 
given a whole-genome resequencing data set for a human with specific diseases, 
typically around 3 million SNPs and around half million insertions/deletions 
will be identified. Given this massive amounts of data (and candidate disease-
causing variants), it is necessary to have a fast algorithm that scans the data 
and identify a prioritized subset of variants that are most likely functional 
for follow-up Sanger sequencing studies and functional assays.

Currently, these various types of functional annotations produced by ANNOVAR can 
be (1) gene-based annotations (the default behavior), such as exonic variants, 
intronic variants, intergenic variants, downstream variants, UTR variants, 
splicing site variants, stc. For exonic variants, ANNOVAR will try to predict 
whether each of the variants is non-synonymous SNV, synonymous SNV, 
frameshifting change, nonframeshifting change. (2) region-based annotation, to 
identify whether a given variant overlaps with a specific type of genomic 
region, for example, predicted transcription factor binding site or predicted 
microRNAs.(3) filter-based annotation, to filter a list of variants so that only 
those not observed in variation databases (such as 1000 Genomes Project and 
dbSNP) are printed out.

Currently, I am expanding the functionality of ANNOVAR on (1) Fusion gene 
detection from large deletions, where a deletion joins the reading frame of two 
genes (same orientation of transcription) together to create a new gene. (2) 
Assignment of functional importance score to each observed mutation in the 
genome. This will be extremely important for the development of association 
tests for rare variants, and for prioritization of variants in downstream 
functional studies after a successful genome-wide association studies (GWAS).

=over 8

=item * B<variant file format>

A sample variant file contains one variant per line, with the fields being chr, 
start, end, reference allele, observed allele, other information. The other 
information can be anything (for example, it may contain sample identifiers for 
the corresponding variant.) An example is shown below:

	16      49303427        49303427        C       T       rs2066844       R702W (NOD2)
	16      49314041        49314041        G       C       rs2066845       G908R (NOD2)
	16      49321279        49321279        -       C       rs2066847       c.3016_3017insC (NOD2)
	16      49290897        49290897        C       T       rs9999999       intronic (NOD2)
	16      49288500        49288500        A       T       rs8888888       intergenic (NOD2)
	16      49288552        49288552        T       -       rs7777777       UTR5 (NOD2)
	18      56190256        56190256        C       T       rs2229616       V103I (MC4R)

=item * B<database file format: UCSC Genome Browser annotation database>

Most but not all of the gene annotation databases are directly downloaded from 
UCSC Genome Browser, so the file format is identical to what was used by the 
genome browser. The users can check Table Browser (for example, human hg18 table 
browser is at http://www.genome.ucsc.edu/cgi-bin/hgTables?org=Human&db=hg18) to 
see what fields are available in the annotation file. Note that even for the 
same species (such as humans), the file format might be different between 
different genome builds (such as between hg16, hg17 and hg18). ANNOVAR will try 
to be smart about guessing file format, based on the combination of the --
buildver argument and the number of columns in the input file. In general, the 
database file format should not be something that users need to worry about.

=item * B<database file format: GFF3 format for gene-based annotations)>

As of June 2010, ANNOVAR cannot perform gene-based annotations using GFF3 input 
files, and any annotations on GFF3 is region-based. However, this is expected to 
be changed in the future.

=item * B<database file format: GFF3 format for region-based annotations)>

Currently, region-based annotations can support the Generic Feature Format 
version 3 (GFF3) formatted files. The GFF3 has become the de facto golden 
standards for many model organism databases, such that many users may want to 
take a custom annotation database and run ANNOVAR on them, and it would be the 
most convenient if the custom file is made with GFF3 format. 

=item * B<database file format: generic format for filter-based annotations)>

The 'generic' format is designed for filter-based annotation that looks for 
exact variants. The format is almost identical to the ANNOVAR input format, with 
chr, start, end, reference allele, observed allele and scores (higher scores are 
regarded as better).

=item * B<database file format: VCF format for filter-based annotations)>

The 1000 Genomes Project now provide their variant annotations in VCF format, so 
I implemented the functionality to directly interrogate VCF files. A VCF file 
may contain summary information for variants (for example, this variant has MAF 
of 5% in this population), or it may contain the actual variant calls for each 
individual in a specific population. As of March 2010, the files from 1000G website 
only contains the first type of information (that is, alleles and their 
frequencies in population). For the purpose of simplicity, ANNOVAR only 
interrogates the first type of information.

=item * B<database file format: avsift for filter-based annotations)>

avsift refers to a file that ANNOVAR developers compiled for fast annotation of 
SIFT scores for non-synonymous variants in the human genome. It conforms to the 
generic format described above. However, users can directly specify '--dbtype 
avsift' in command line to perform avsift annotations, making it more convenient 
for users. Alternatively, users can use '--dbtype generic -genericdbfile 
hg18_avsift.txt' for the annotation, and the effects are usually the same.

=item * B<sequence file format>

ANNOVAR can directly examine FASTA-formatted sequence files. For mRNA sequences, 
the name of the sequences are the mRNA identifier. For genomic sequences, the 
name of the sequences in the files are usually chr1, chr2, chr3, etc, so that 
ANNOVAR knows which sequence corresponds to which chromosome. Unfortunately, 
UCSC uses things like chr6_random to annotate un-assembled sequences, as opposed 
to using the actual contig identifiers. This causes some issues (depending on 
how reads alignment algorithms works), but in general should not be something 
that user need to worry about. If the users absolutely care about the exact 
contigs rather than chr*_random, then they will need to re-align the short reads 
at chr*_random to a different FASTA file that contains the contigs, and then 
execute ANNOVAR on the newly identified variants.

=item * B<invalid input>

If the query file contains input lines with invalid format, ANNOVAR will skip 
such line and continue with the annotation on next lines. These invalid input 
lines will be written to a file with suffix invalid_input. Users should manually 
examine this file and identify sources of error.

=back

-----------------------------------------

ANNOVAR is free for academic, personal and non-profit use.

For questions or comments, please contact kai@openbioinformatics.org.

=cut
