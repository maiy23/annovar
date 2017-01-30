#!/usr/bin/perl

#THIS IS A TEST EDIT
use strict;
use FindBin;
=head
This script is a wrapper for fdi
=cut
use Data::Dumper;
 umask(0000);
$ENV{'ORACLE_HOME'}='/opt/nasapps/production/oracle/product/11.2.0/client';
my $usage=qq{
	$0 -i <dbname> <comma separated list of variants> <OPTIONAL Organism>
		Note: Variants should be in the form: 
			chrX:XXXX:XXXX:REF:VAR 
		e.g.
			chr1:10611:10611:C:G
		Note: Organism should be in the form hg19 or mm9
		
};
my $proc_dir='/bioinfoC/AVA/FDI/';
if ($#ARGV<1){print "$usage\n";exit;}
my $fdi_dbname=$ARGV[0];
my $opt_o= ($#ARGV>1)?$ARGV[2]: 'hg19';
my $orgdb="humandb";
if ($opt_o!~/hg/i){$orgdb="mousedb";}
my $SEQIDX_DIR ="/SeqIdx/annovardb/$orgdb";
my $cwd=`pwd`;chomp $cwd;
my $bin=`dirname $0`;chomp $bin;
die "[ERROR] Your database directory $SEQIDX_DIR does not exist or is not a directory\n"  if ( (defined $SEQIDX_DIR) && (!-e $SEQIDX_DIR || !-d $SEQIDX_DIR) );
my @files;
#### make the file for ANNOVAR input
my $input_filename="testmirna";#`date "+%Y%m%d%S"`;chomp $input_filename;
$proc_dir.="$input_filename";
system ("mkdir -p $proc_dir\n");
chdir ("$proc_dir");
open (FILE,">$input_filename") or die "Cannot open $input_filename\n";
foreach my $var (split(",",$ARGV[1])){
	$var=~s/:/\t/g;
	next if ($var eq '');
	print FILE "$var\n";
}
close FILE;
my $fetch_script="$bin/../annofetch";
my ($addon,$cmd2run);
my %dummy_hash;
#print "\n================\n$fdi_dbname\n================\n";#for testing
my $outname;my $notgrid=0;
if ($fdi_dbname=~/geneanno/i){die 'need upload schema';
	if ($fdi_dbname=~/ens/){
		$addon="-E";
	}
	$cmd2run="perl $bin/annovar_qsub_wrpr.pl -i $input_filename $addon";
}elsif ($fdi_dbname=~/miRNA_abcc/i){
	#first put it into the needed format for Natalia's script
	open (BAT,">mirna.bat") or die "cannot open bat file\n";
	print BAT "#PBS /bin/bash\n#PBS_O_WORKDIR=$proc_dir\n#PBS -j oe -e $proc_dir/annovar_stderrout\n";	
	print BAT "cd $proc_dir/;perl $bin/annotate_variation_ABCC2.pl --regionanno --dbtype gff3 --gff3dbfile hsa.gff3  --keepline --noheader $input_filename /bioinfoC/hue/mirna/  \n";
	close BAT;
	system ("chmod +x mirna.bat");
#	runAndMonitor("qsub mirna.bat",__LINE__);
	#rewrite output
	open (FILE,"<$input_filename.$opt_o\_gff3 ") or die "Cannot open file $input_filename.$opt_o\_gff3 \n";
	open (OUT,">tmp.mirna.input" ) or die "cannot open tmp.mirna.input\n";
	my @names;
	while (my $line=<FILE>){
		chomp $line;
		next if ($line=~/^[\t\s]*$/);
		next if ($line=~/^#/);
		my ($name,@rest)=split("\t",$line);
		if ($name=~/Name=(\S+)/){
			$name=$1;
			@names=split(",",$name);
		}else{
			$dummy_hash{$rest[1]}{'id'}=join(":",@rest);
			$dummy_hash{$rest[1]}{'impact'}="-";
			next;#do not submit if unknown
		}
		$dummy_hash{$rest[1]}{'id'}=join(":",@rest);
		
		foreach my $name (@names){
			$dummy_hash{$rest[1]}{$name}=0;#hsa name and chrpos as key
			print OUT join ("\t",$name,@rest)."\n";
		}
	}
	close OUT;
	$outname="Impact_Table";
	if (!-z "tmp.mirna.input"){
#		$cmd2run="perl /bioinfoC/hue/mirna/mirna_impact_all.pl $proc_dir tmp.mirna.input";
	}else{
#		$cmd2run ="touch $proc_dir/$outname";$notgrid=1;
	}
}elsif ($fdi_dbname=~/flank/i){
	print STDERR "Fetch script did not work..skipping ($fetch_script)" and exit if (! -e $fetch_script);
	$cmd2run="$fetch_script $input_filename $input_filename.flanks";
}else{
	open (DB,">$input_filename.database") or die "Cannot open $input_filename.database\n";
	print DB "$SEQIDX_DIR/$opt_o\_$fdi_dbname.txt";
	close DB;
	$cmd2run="perl $bin/annovar_qsub_wrpr.pl -f $input_filename.database -i $input_filename -N";#capture the grid submission id
}
#runAndMonitor("$cmd2run",__LINE__);
my $addon;
if ($fdi_dbname=~/mirna_abcc/i){
	#format of Impact_Table:
	#hsa-mir-604     ->      hsa-mir-604_-_0:29833998;T/C    pseudo_miRNA(67.8%)->not_miRNA()        affected
	open (FILE,"<$outname") or die "Canot open!\n";
	while (<FILE>){
		next if ($_=~/^#/);
		my ($mirna,$arrow,$info,$compare,$impact)=split("\t",$_);chomp $impact;
		my $pos;
		if ($info=~/\_[\+\-]\_\d+:(\d+);\S+/){
			$pos= $1;
		}else{
#			print "skipping $info\n";
			next;
		}
		if (exists $dummy_hash{$pos}{$mirna}){
			$dummy_hash{$pos}{'impact'}.="$mirna:IMPACT=$compare:$impact,";
		}else{
			die "doesnot exist $pos $mirna\n";
		}
	}
	close FILE;
	open (OUT, "|tee $input_filename.annovar_wrpr.output") or die "Cannot open $input_filename.annovar_wrpr.output for writing in $0\n";
	$outname="$input_filename.annovar_wrpr.output";
	foreach my $pos (sort keys %dummy_hash){
		print OUT "$dummy_hash{$pos}{'id'}\t";
		if ($dummy_hash{$pos}{'impact'}){
			$dummy_hash{$pos}{'impact'}=~s/^0//;
			chop $dummy_hash{$pos}{'impact'} if ($dummy_hash{$pos}{'impact'}=~/,$/);
			print OUT "$dummy_hash{$pos}{'impact'}\n";
		}else{
			print OUT "-\n";
		}
	}
	close OUT;
	$addon="$outname.out 0 $fdi_dbname";#default name;
}else{
	$outname="$input_filename.annovar_wrpr.output";
	if (-e $outname){
		open (FILE,"<$outname") or die "Canot open!\n";
		while (<FILE>){
			next if ($_=~/^#/);
			my ($value,@others)=split("\t",$_);chomp $others[$#others];
			print STDOUT join (":",@others) . "\t$value\n";
		}
		close FILE;
	}
}
exit;#for now we are not loading the database
if (-z "$input_filename.annovar_wrpr.output"){
	exit;
}
system ("perl $bin/parse4db.pl  $input_filename.annovar_wrpr.output $addon\n");
my $ctl=qq(OPTIONS (DIRECT=TRUE, BINDSIZE=2000000)
UNRECOVERABLE
LOAD DATA
INFILE '$input_filename.annovar_wrpr.output.out'
APPEND
INTO TABLE hg19_impacts 
FIELDS TERMINATED BY X'09' 
(chr, var_start, var_end, ref_allele char(200) "substr(:ref_allele,0,200)", var_allele, impact_name, impact_value char(2500) "substr(:impact_value,0,2500)")
);
open (OUT,">impacts.ctl") or die "cannot open file";
print OUT "$ctl\n";
close OUT;
system ("$ENV{ORACLE_HOME}/bin/sqlldr vuonghm/hue0404\@abcc11d control=impacts.ctl\n");
exit();
sub runAndMonitor{
	my $cmd=shift;
	my $line=shift;
	my $grid_id=`$cmd`;
	if ($grid_id=~/(\d{6,})/){
		my $targets="qstat | grep $1 ";
		my $timeout=`date "+%d%H%M%s"`;chomp $timeout;$timeout+=6000;#1hour
		my $done=1;
		until (`$targets` eq '' || !$done){
			sleep (15);
			my $curr_time=`date "+%d%H%M%S"`;chomp $curr_time;
			if ( $curr_time>$timeout){$done=0;}
		
		};
		return if ($done);
		die "Could not finish the run ($targets) for $cmd at $line\n";
	}elsif (!$notgrid){
		#reset each time;  global variable; set at command specification
		$notgrid=0;
		return 1;
	}else{
		die "Could not run $cmd and ($grid_id)\n";
	}
}
