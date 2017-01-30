#!/usr/bin/perl
use DBI;
use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use reseq_modules2;
use vars qw( $opt_v $opt_P $opt_d $opt_t $opt_o $opt_h $opt_H $opt_l $opt_e $opt_p);
getopts("v:d:P:p:o:l:thHe");
my $SIFT_HOME='/bioinfoA/sift';
my $usage=qq{
	$0
		-P <filename with positions> TAB DELIMITED 
			ONE VAR PER LINE
			CHR->POSITION(spacebased)->Orientation->Allele
			13->31804729->1->A/C
		-p <ANNOVAR input>
			ONE VAR PER LINE 
			CHR->POS1->POS2->REF1->ALLELE2
	OPTIONAL:
		-v <database version> DEFAULT:37
		-H <print Header>
		-e <no fail> 
			useful when running from another script and especially long lists
			prints error message to LOG file
		-o <output filename> DEFAULT : prints to STDOUT
		-d <delimiter>	DEFAULT: tab
			changes the delimiter from tab to user input
		-l <log file>
		-h <prints this help message>
}; 
if ($opt_h){
	print $usage;
	exit;
}
my $start_time=timings(0);
my $user="abccruser";#sqltest3 user'ABCC_GB_ROUser';
my $pass="abccrpwd";#sqltest3 pass"abcC";
die "Please specify a positions file($opt_P).\n$usage\n" if (!defined $opt_P && !defined $opt_p);
my $chr; #for example
my $ver=37;
my $log_fn;
my $inputfile=(defined $opt_P)?$opt_P:$opt_p;
my $cwd=`dirname $inputfile`;chomp $cwd;

if (defined $opt_l){
	$log_fn=$opt_l;
}else{
	$log_fn="$cwd/sift_mysql_wp.log";
}
open (LOG, ">$log_fn") || die "Cannot open $log_fn\n";
print LOG "[INFO] Started at " . `date`;
if (defined $opt_v){
	$ver=$opt_v;
	die "ERROR in $0 and the version\n" if ($ver!~/\d+/);
}
print LOG "Using v$ver databases...\n";
my $delimiter;
if ($opt_d){
	$delimiter=$opt_d;
}else{
	$delimiter='\t';
}

my @pos;my %query; my %indels;
if ($opt_P){
	die "File ($opt_P) was not defined or does not exist\n\n$usage\n" if (!defined $opt_P || !-e $opt_P);
	open (FILE,"<$opt_P") || die "Cannot open your inputfile ($opt_P)\n";
	while (<FILE>){
		next if ($_ eq '');
		chomp;
		my @info=split("$delimiter",$_);
		if ($info[0]=~/chr([XY\d]+)/){
			$info[0]=$1;
		}
		if ($_=~/(ins|del)/){
			$indels{"$info[0]\_$info[1]"}{'info'}=join ",",@info;
			push(@{$query{"$info[0]\_$info[1]"}},join (",",@info))
		}elsif (exists $query{"$info[0]\_$info[1]"}){
			push(@{$query{"$info[0]\_$info[1]"}},join (",",@info)) ;
		}else{
		 	push(@{$query{"$info[0]\_$info[1]"}},join (",",@info))
		}
	}
	close FILE;
}elsif ($opt_p){
	die "File ($opt_p) was not defined or does not exist\n\n$usage\n" if (!defined $opt_p || !-e $opt_p);
	open (FILE,"<$opt_p") || die "Cannot open your inputfile ($opt_p)\n";
	while (<FILE>){
		next if ($_ eq '');
		chomp;
		my @info=split("$delimiter",$_);
		if ($info[0]=~/chr([XY\d]+)/){
			$info[0]=$1;
		}
		if ($info[3] eq '-' || $info[4] eq '-'){
			$indels{"$info[0]\_$info[1]"}{'info'}=join ",",@info;
		}elsif (exists $query{"$info[0]\_$info[1]"}){
			push(@{$query{"$info[0]\_$info[1]"}},join (",",@info)) ;
		}else{
			push(@{$query{"$info[0]\_$info[1]"}},join (",",@info))
		}
	}
	close FILE;
}else{
	die "you must specify the correct input file and format\n";
}
print LOG "[INFO] Done reading input file ($opt_P)\n";
print LOG  Dumper (\%query,"LNNBR".__LINE__);
my $interval=timings($start_time);
print LOG "[INFO] Reading in input and printing Dumper $interval\n";
my $table_chr="";
my %bins;#this is to bin all of the queries so that the tables can queried
my $dbh= DBI->connect("dbi:mysql:SIFT_Human_$ver:sqldb2.abcc","$user","$pass") || die "Could not connect to the database\n";
#find the associated table that will eventually query
my $sql="select chrbin_chrom,chrbin_start,chrbin_stop from chrom_bins where chrbin_chrom=? and chrbin_start<=? and chrbin_stop>=?";
my $sth=$dbh->prepare ($sql);
print LOG "preparing $sql statement";
foreach my $qry (keys %query){
	my ($chr,$pos)=split ("_",$qry);my ($start,$stop);
	if ($pos=~/(\d+)\-(\d+)/){
		if ($1<$2){
			$start=$1;$stop=$2;
		}else{
			$stop=$1;$start=$2;
		}
	}else{
		$start=$stop=$pos;
	}
	print LOG "Running $chr,$start,$stop\n";
	$sth->execute ($chr,$start,$stop);
	while (my @results=$sth->fetchrow_array()){
		my $bin=join ("_",@results);
		print LOG "[INFO] BINNBR:$bin\n";
		push (@{$bins{$bin}},$query{$qry});
	}
	if ($sth->rows==0){
		print LOG "\n\n\n\nERROR\n";
		die "ERROR";
	}
}
print LOG "[INFO] Done finding bins\n";
undef (%query);
#
##(CHR TEXT NOT NULL , COORD1 NUMERIC NOT NULL , COORD2 NUMERIC NOT NULL , ORN TEXT, RSID TEXT, ENSG TEXT, 
##ENST TEXT, ENSP TEXT, REGION TEXT, SNP TEXT, NT1 CHAR, NT2 CHAR, NTPOS1 NUMERIC, NTPOS2 NUMERIC, CODON1 TEXT, 
##CODON2 TEXT, AA1 CHAR, AA2 CHAR, AAPOS1 NUMERIC, AAPOS2 NUMERIC, CDS NUMERIC, AA1_VALID INTEGER, 
##ENST_VALID INTEGER, SCORE NUMERIC, MEDIAN NUMERIC,  INTEGER, PRIMARY KEY (CHR, COORD2, ORN, ENST,NT2,AA2))
my $output="";
if (defined $opt_o){
	open (OUT,">$opt_o") || die "Cannot open $opt_o for writing\n";
}else{
	open (OUT,">&STDOUT") || die "cannot open STDOUT!!\n";
}
my %predict;#print Dumper (\%indels);
my $indel_sth;my $count=0;
if ($opt_H){print "CHR,COORD,VAR_OBS,ENST,NTPOS1,NTPOS2,AA1,AA2,AAPOS1,AAPOS2,SCORE,Median,REFBC\n";}
print LOG "[INFO] Ready to print\n";
$interval=timings($start_time);
print LOG "[INFO] $interval\n";
print LOG Dumper (\%bins,"printed binshash");
foreach my $table_chr (keys %bins){
	my ($chr,@range)=split("_",$table_chr);$chr="chr$chr";
	my $query="select coord2,nt2,enst,ntpos1,ntpos2,aa1,aa2,aapos1,aapos2,score,median,nt1 from chr$table_chr where coord2=? and (nt2=? || nt2=?)";
#	my $query="select * from chr$table_chr where coord2=? and (nt2=? || nt2=?)";
	my $indel_qry="select coord2,nt2,enst,ntpos1,ntpos2,aa1,aa2,aapos1,aapos2,score,median,nt1 from chr$table_chr where (? between coord1 and coord2 || ? between coord1 and coord2)";
	$sth=$dbh->prepare($query);
	$indel_sth=$dbh->prepare ($indel_qry);
	print LOG "[INFO] working on $table_chr ".$#{$bins{$table_chr}}." queries\n";
	foreach my $xinfo_arr (@{$bins{$table_chr}}){#ARRAY OF ARRAYS of coordinates [ [1,1234,1,A/C], [ [ '2,1234,-1,A/G','2,1234,-1,A/C'] };
		foreach my $query(@{$xinfo_arr}){#array of one variant with chr,pos,orn,alleles
			my ($chr,$pp,$orn,$nt)=split(",",$query);
			next if ($pp eq '');$count++;
			if ($count%100==0){print LOG "\tQuery:$count\t$query\n";}
			if ($nt!~/(ins|del)/){
				$nt=uc($nt);
				my ($nt1,$nt2)=split("/",$nt);
				if ($orn=~/\-/){
					$nt1=~tr/[ATCG]/[TAGC]/;
					$nt2=~tr/[ATCG]/[TAGC]/;
				}
				$sth->execute($pp,$nt1,$nt2) or go_die ("Couldn't execute statement: at LN145 " . $sth->errstr,0);
				while (my @rows = $sth->fetchrow_array()){
					if ($orn=~/\-/){
						$rows[1]=~tr/[ATCGatcg]/[TAGCtagc]/;
						$rows[11]=~tr/[ATCGatcg]/[TAGCtagc]/;
					}
					$output.= join (",",$chr,$pp,@rows[1..$#rows],"\n");
				}
				if ($sth->rows == 0) {
		          $output.="$chr,$pp,$nt,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA\n";
		       }
			}else{
				if ($pp=~/(\d+)\-(\d+)/){
					$indel_sth->execute ($1,$2);
				}else{
					$indel_sth->execute ($pp,$pp);
				}
				while (my @rows = $indel_sth->fetchrow_array()){
					my @subset=@rows[2..8];
					 if (!exists $indels{"$chr\_$pp"}{$subset[0]}){
					 	$indels{"$chr\_$pp"}{$subset[0]}=join (",",@subset);
					 	$indels{"$chr\_$pp"}{'ref'}=$rows[11];
					 }
				}
				if ($indel_sth->rows == 0) {
					foreach my $element (keys %{$indels{"$chr\_$pp"}}){
						undef $indels{"$chr\_$pp"}{$element};
					}
					delete $indels{"$chr\_$pp"};
		           $output.="$chr,$pp,$nt,NA,NA,NA,NA,NA,NA,NA,NA,NA,\n";
		       }
			}
			if ($count%2400 ==0){
				print LOG "\t\tDUMPING $count\n";
				print OUT "$output";$output="";
			}
		}
	}
	delete ($bins{$table_chr});
	print OUT "$output";
	$output="";
}
$interval=timings ($start_time);
print LOG "[INFO] $interval\n";
my %mapper;

if (keys %indels ==0){
	print LOG "[INFO] Completed!\n";exit;
}
my $sift_tmp_fn;
if ($opt_P=~/sift.inp.(\d+)/){#specified pid
	$sift_tmp_fn="$cwd/tmp.indel.$1";
}else{
	#we do this in case another program is also running this script in the same directory 
	#ie. parallel jobs
	my $ct=1;
	until (!-e "$cwd/sift.tmp.indel.$ct"){
		$ct++;
	}
	$sift_tmp_fn="$cwd/sift.tmp.indel.$ct";
}
open (INDELS,">$sift_tmp_fn") or die "Cannot open ` . $!\n";
foreach my $indel (keys %indels){
	if (!exists $indels{$indel}{'info'}){next;}
	my ($chr,$pos,$orn,$type)=split(",",$indels{$indel}{'info'});
	if (($chr=~/(\d+)\.(\d+)/) && $chr!~/\-/){
		$chr=$1;
	}
	my $range=$pos;my $siftinpact;$range=~s/-/,/g;
	if ($type=~/ins(.*)/){
		print INDELS "$chr,$range,$orn,$1\n";
	}elsif ($type=~/del/){
		print INDELS "$chr,$range,$orn,-/\n";
	}else{
		print "did not catch this($indel)\n";
	}
	$mapper{"$indel"}="$type";
} 	
if (!-e "$cwd/sift_dir"){
	mkdir ("$cwd/sift_dir");
}
print LOG "[INFO] About to run:/bioinfoA/sift/sift_mysql/bin/SIFT_exome_indels.pl -c /bioinfoA/sift/db/Coding_info_$ver -i $sift_tmp_fn -d SIFT_Human_$ver -o $cwd/sift_dir 2>&1|head -n1\n";
my $pid=`/bioinfoA/sift/sift_mysql/bin/SIFT_exome_indels.pl -c /bioinfoA/sift/db/Coding_info_$ver -i $sift_tmp_fn -d SIFT_Human_$ver -o $cwd/sift_dir 2>&1|head -n1`;
chomp $pid;
if ($pid=~/Your job id is (\d+) and is currently/){
	$pid=$1;
}else{
	go_die( "Cannot parse PID ($pid) after running: at LN224\n".
	"/bioinfoA/sift/sift_mysql/bin/SIFT_exome_indels.pl -c /bioinfoA/sift/db/Coding_info_$ver -i $cwd/$sift_tmp_fn -d SIFT_Human_$ver -o $cwd/sift_dir |head -n1\n");
}
my $rtn="";my $hasrpt="";my $wait_ctr=0;
my $max=keys %indels;my $max_counter;

#timeout counter
if ($max%100<=1){
	$max_counter=100;
}else{
	$max_counter=$max*100;
}
print LOG "[INFO] Waiting for sift_indels for $pid\n";
until ((-e "$cwd/sift_dir/$pid/$pid.classification_file")||($wait_ctr>$max_counter)){
	sleep(20);$wait_ctr++;
}
if ($wait_ctr>$max_counter){print LOG "waiting to long for sift\n";go_die( "ERROR with SIFT_INDEL at LN240 sub $0\n");}
if (-e "$cwd/sift_dir/$pid/$pid.classification_file" && !-z "$cwd/sift_dir/$pid/$pid.classification_file") {
	open (CF, "$cwd/sift_dir/$pid/$pid.classification_file") || go_die ("LN 242:cannot open $cwd/sift_dir/$pid/$pid.classification_file\n");
	print LOG "opening $cwd/sift_dir/$pid/$pid.classification_file\n";my $flag=0;my $lastpos="";
	while (<CF>){
		my @arr=split("\t",$_);
		my ($c,$pos1,$pos2,$orn,$value)=split (",",$arr[0]);
		my ($indel_pos,$indel);
		if (($pos1==$pos2) && ($value eq "/")){
			$indel_pos=$pos1;
			$indel="het_indel(<3)";
		}elsif ( ($pos1==$pos2) && ($value=~/^([ATCG]*)$/) ){
			my $len=length($1);$indel="ins$1";
			$indel_pos="$pos1.1-$pos1.$len";
		}elsif ( $pos1!=$pos2 ) {
#				die "not tested yet in ($pos1 and $pos2 )$0 sift_indel subat LN". __LINE__."\n";
			$indel_pos="$pos1-$pos2";
			if ($arr[6]=~/\-([ATCG]*)\-/i){
				$indel=uc($1);
				if ($orn eq '-'){$indel="del".revcmpl($indel);}else{$indel="del$indel";}
#print "Found $indel AT $indel_pos --> $pid\n";my $wait=<STDIN>;
			}else{
				die "can't parse this ".join (",",@arr)."\n";
			}
		}elsif ($pos1==$pos2){
			$indel_pos="$pos1-$pos2";
		}else{
			die "dont know how to parse $c\_$pos1-$pos2\n";
		}
		if (! exists $indels{"$c\_$indel_pos"}{$arr[3]}){go_die ( "$c\_$indel_pos and $arr[3] does not exist LN 268\n");next;}
		if ($arr[8]=~/intron/i){
			$rtn="NA,NA,NA,NA,NA,NA,NA,NA,NA";
			print OUT "$c,$indel_pos,".$mapper{"$c\_$indel_pos"}.",$arr[3],NA,NA,NA,NA,NA,NA,NA,NA,".$indels{"$c\_$indel_pos"}{'ref'}."\n";
		}else{
			my $perc;
			$rtn="";
			if ($arr[13]=~/no/i){$rtn.="Not NMD";}elsif ($arr[13]=~/yes/i){$rtn.="NMD";
			}else{
				#http://sift.jcvi.org/www/chr_coords_example_indels.html
				#There is no NMD when:
					#1) the resulting premature termination codon is in the last exon
					#-or-
					#2) the resulting premature termintion codon is in the last 50 nucleotides in the second to last exon 
				$flag=1;
				if ($arr[11]=~/(\d+)\%/){
					$perc=$1;
				}
				if ($perc<=80){
					$rtn.="Pred NMD";
				}else{
					$rtn.="Pred Not NMD";
				}
			}
			if ($arr[15] ne ""){
				$hasrpt=$arr[15];
			}
			$rtn.= " ($arr[11])";
			print OUT "$c,$indel_pos,".$mapper{"$c\_$indel_pos"}.",".join (",",$indels{"$c\_$indel_pos"}{$arr[3]}).",$rtn,,".$indels{"$c\_$indel_pos"}{'ref'}."\n";#print Dumper (\%indels);die;
		}
	}
	close CF;
	
}elsif (-e "$cwd/sift_dir/$pid/$pid.classification_file" && -z "$cwd/sift_dir/$pid/$pid.classification_file") {
	foreach my $indel_pos(split ("\n",`cut -f2 -d, $sift_tmp_fn`)){
		print "working on sift_indel $indel_pos\n";
#			$$hash_ref{$indel_pos}{$indel}=0;
#		$impact{$indel_pos}{'sift'}="Mutation is in an mrna (intron)";
	}
}else{
	go_die ( "Could not find pid $cwd/sift_dir/$pid/$pid.classification_file\n",0);
}
$interval=timings ($start_time);
print LOG "[INFO] $interval\n";
print LOG "[INFO] Completed\n";


exit;
sub go_die{
	my $msg=shift;
	my $dumper=shift;
	if (ref($dumper)=~/HASHREF/i){
#		print Dumper LOG (\%{$dumper});
	}
	print LOG "[DATAERR] $msg\n";
	return;
}