#!/usr/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin/Parallel-ForkManager-0.7.9/lib";
use Parallel::ForkManager;
use List::MoreUtils 'pairwise';#to add arrays of same size together by element
=head
This script is a wrapper for annovar and its submission to the grid
outputs all the necessary files for grid 

You can specify a config file so that you can control your annovar arguments,
but mostly it works as is.
v2.0 Add batching for parallelization
	ForkManager pm used for parallelization of the batch mode (-b)
v2.1 Corrected for user specified -g and running on tork
	fixed call to annotate_variation.pl so it uses $bin and not hard-coded full path.
	automatically batches data (-b is no longer valid)
v3.0 error handling for when batch(es) fail.  need to alert user of internal error and that their data is not complete.
	--diagnose fail type and rerun? (planning)
=cut
if (join("",@ARGV)=~/\-W/  || ($ENV{'HOST'}=~/tork/i && $#ARGV>-1 && join("",@ARGV)!~/\-g/)){
	#do not run CPU intensive jobs on tork
	my $dir=`pwd`;chomp $dir;$dir=`basename $dir`;chomp $dir;
	open (MASTER, ">PBS.MASTER") or die "cannot open PBS.MASTER for writing\n";
	print MASTER "#!/bin/sh\n#PBS -j oe -l arch=linux-x86_64,pvmem=2GB,mem=2GB,ncpus=2 -m ea\n";
	print MASTER "cd ".`pwd`;
	my $cmd= "perl $0 ".join (" ",@ARGV)." -g -T 2\n"; 
	$cmd=~s/\-W//;
	print MASTER "$cmd";
	close MASTER;
	exec ("qsub -q grande -N PBS_$dir PBS.MASTER\n");
	exit;
}else{
	print "$0 ".join (" ",@ARGV)."\n";
}
use Getopt::Std;
use Data::Dumper;
use FindBin;
use Tie::File;
use lib "$FindBin::Bin/";
use ABCC;
use vars qw($opt_W $opt_I $opt_E $opt_R $opt_A $opt_G $opt_w $opt_z $opt_F $opt_x $opt_m $opt_H $opt_f $opt_a $opt_o $opt_i $opt_d $opt_T  $opt_h $opt_t $opt_c $opt_g $opt_n $opt_l $opt_b $opt_B $opt_k $opt_r $opt_s $opt_N $opt_v $opt_S $opt_p);
 getopts("f:o:c:d:i:T:t:l:B:k:p:m:F:G:A:IrvzsbcaehxgnHINSwRWE");
 umask(0000);
my $usage=qq{
	$0 -i <input filename>
	[ONE OF THE FOLLOWING MUST BE SPECIFIED]
		-f <filename with databases to use>
		-a <if specified, then it runs through all the databases for that organism in the directory>
	[OPTIONAL]
		-g <if specified, do not run on grid>
		-o <organism>  hg19
		-d <full path to database directory> 
			DEFAULT: /SeqIdx/annovardb
		-t <txt|vcf>
		-c <if specified then comments will not be included>
		-T <length of time before timeout> DEFAULT 1hr
			in units of hours
		-n <no run, just aggregate data>
				-use the -k <num> option to specify the number of files that you wish your data to be split
		-l <if aggregate was interrupted, use this to specify the line in the input file>
			automatically specifies -n option
		-b <if specified, then all will run in parallel blocks>
			ie. if you have a file of 1,000,000, this script will split the outputs into blocks of 10,000 (i.e. 10 files) so you can run them in parallel.
				It will also aggregate each of those files into one file at the end.
			You may specify -B option to specify the size of your block.  DEFAULT=1,000,000
		-r <if specified, then it does NOT delete intermediate files>
		-p <output filename> 
			#default: annovar_wrpr.output 
		-x <do not use strict>
			If possible, AVIA uses filter annot option to use the exact alleles to annotate, not allowing for any base
			(filter instead of regionanno)
			if specified, it will turn OFF this functionality and be more lenient for annotation
		-m <column for name for indexing>  ##ADVANCED OPTION for AVIA pipeline, 0based
			Only used if indexing occurs on the input database file.  
		-F <file with filtering order and instructions on which files to keep>
			input file looks like (one per line)
				(1000g|hapmap|snp135|pp2|siftv63)\t(filtered|dropped|dropcutoff|keepcutoff)\tcutoff value (if cutoff)
		-z <if specified, then also runs query to get flanking sequences>
		-G <specify a gene list in which to filter your output>
		-A <allele frequency file>
		-R <rename headers and adds the filename to the beginning column>
			#renames the headers of the output file to a more reader friendly format
			#also adds the filename as a column at the beginning of the file e.g. name_header_column=1 in config file
		-I <runs NMD Indel script>
		-W <is a web input>   #reduces the number of nodes allowable for use 
	[CONDITIONAL]
		-w if specified, then it will generate web files;
		-k <number of files to split your input files>
			this is only used for -b option (batches).  Specifies the number of batches.
		-s <do not split files>
		-N <do not run Annotation>
		-v <if specified then the input file only has chr,start and stop positions, exonic function will not be run
		-S <if specified, sorts the input so that it *may* run faster if using indexed files>
		
		
};
if ($opt_h || !defined $opt_i){print "$usage\n";exit;}
$opt_o ||= 'hg19';
$opt_g ||= '0';
$opt_B ||= 2_000_000;#in tests this seems to have the same runtime for batched vs non-batched, if not faster than the batched 2/06/12
$opt_T = (defined $opt_T)? $opt_T*10000:60000;
if (defined $opt_b){
	$opt_T=60000;
}
if (!$opt_c){$opt_c=" --comment ";}else{ $opt_c='';}
$opt_t ||= "txt";
my $out =(defined ($opt_p))?"$opt_p.annovar_wrpr":"$opt_i.annovar_wrpr";
my %rank;#this is your rank of most important impact to least important
my $orgdb="humandb";
if ($opt_o!~/hg/i){$orgdb="mousedb";}
$opt_d ||="/SeqIdx/annovardb/$orgdb";
my $cwd=`pwd`;chomp $cwd;
if (defined $opt_l){$opt_n=1;print STDERR "[INFO]  You are specifying line $opt_l...using -n as well..\n";}
my $bin=`dirname $0`;chomp $bin;
my $annovar="perl $bin/annotate_variation_ABCC.pl";
die "[ERROR] Your database directory $opt_d does not exist or is not a directory\n"  if ( (defined $opt_d) && (!-e $opt_d || !-d $opt_d) );
die "[ERROR] Your -k $opt_k entry was invalid\n$usage" if (defined $opt_k && $opt_k!~/^\d+$/);
my @files;
my $MAX_PROCESSES=4;#max number of spans for a parallel process
checkPBSAccess() if (!$opt_n);
my %flg;#keeps track of different flags
my $convertible='FOO';#this is the tag that means that there was a conversion betwen two versions of the genome
################## TESTING ###################
#filtering
my $filter_cmd;
if (defined $opt_F){
	#this is for cascade filtering
	#order of file does not matter
	open (FILTERDB,"<$opt_F") or die "Cannot open filter database\n";
	my $dbdir='/SeqIdx/annovardb/'.$orgdb;
	#fill these in with the order from the web
	my %order=('hapmap_3.3'=>1,
		'1000gALL_sites_2012_04'=>2,
		'avsift'=>3,
		'ljb_pp2'=>4,
		'snp135'=>5
	);
	while (my $line=<FILTERDB>){
		chomp $line;my $threshold;
		next if ($line=~/^\s+$/);
		my ($dbfile,$keep,$cutoff)=split("\t",$line);
		if (exists $order{$dbfile}){
			$order{$dbfile}.=1;print "adding $dbfile from $line\n";
		}else{
			$order{$dbfile}='01';
		}
		if ($dbfile=~/(1000g|snp)/){
			$threshold='--maf_threshold';
		}elsif ($dbfile=~/avsift/){
			$threshold='--sift_threshold';
		}else{
			$threshold='--score_threshold';
		}
		if ($keep=~/(filtered|dropped)/){
			my $input;#we do this to build the query below for input; placeholder for the name based on FOO or BAR
			if ($keep=~/dropped/){
				$input="BAR";
			}else{
				$input="FOO";
			}
			$order{$dbfile}.="$annovar --buildver $opt_o --dbtype $dbfile --filter $input $dbdir\n";
		}elsif($keep=~/keepcutoff/){
			$order{$dbfile}.="$annovar --buildver $opt_o --dbtype $dbfile  --filter $threshold $cutoff --reverse FOO $dbdir\n";
		}elsif($keep=~/dropcutoff/){
			$order{$dbfile}.="$annovar --buildver $opt_o --dbtype $dbfile  --filter $threshold $cutoff FOO $dbdir\n";
		}
		#make sure the next file exists before inserting into next command
		
	}
	my $input=$opt_i;
	foreach my $key (sort {$order{$a} cmp $order{$b} } keys %order){#most efficient way to sort
		$order{$key}=~s/^\d\d//g;$order{$key}=~s/(FOO|BAR)/$input/;
		my $dropped=$1;
		my $dbfile = ($order{$key}=~/\-\-dbtype\s(\S+)\s/)?$1:'';
		next if (!$dbfile);
		$filter_cmd.="$order{$key}";
		$input=($dropped=~/BAR/)?"$input.$opt_o\_$dbfile\_dropped":"$input.$opt_o\_$dbfile\_filtered";
		print "$order{$key}\n";
	}
	my $id=printPBS("filter", '8',$filter_cmd);
	print STDERR "Monitoring filter queries $id\n";
	die "Could not submit $filter_cmd to the grid\n" if (!$id);
	my $notDone=monitor ("qstat | grep $id") ;
	if ($notDone){
		die "Could not filter your stuff\n";
	}
	print STDERR "Continuing with your data\n";
	
	#now transfer your input file to your $opt_i to use
	$opt_i=$input;$flg{'cascade'}++;
}
############################################
#open your input file and read in the databases you wish to search against
if (defined $opt_f){
	open (INFILES,"<$opt_f") or die "Cannot open your input file($opt_f)\nTry specifying the complete path\n";
	@files=<INFILES>;
	close INFILES;
}elsif ($opt_a){
	@files=`ls $opt_d/$opt_o*txt | grep -v refGene.txt`;
}else{
	@files=();#just runs gene anno
}
push (@files, $opt_A) if (defined $opt_A && -e $opt_A);
my $filemem=-s $opt_i;$filemem=int($filemem/100000000);#in GB
#write on PBS job for each...
my $minmem=2+$filemem;#this is the min amount of memory needed to run a small ANNOVAR database, larger files and databases will add more memory
my $jobs_qry='qstat | grep ';
my @gids;#keeps track of the grid ids;
my $use_strict;
$use_strict=(defined $opt_x)?1:0;
#note the output from regulatory is always on a separate line
my $id;my $notDone=1;
#if (-e "$out.output"){ print STDERR "You have 10 seconds to kill this job before your existing $out.output is deleted\n"; sleep (10);system ("rm $out.output headers.output\n");}
if (-e "DONE.$opt_i.anvr.wrpr"){unlink ("DONE.$opt_i.anvr.wrpr");}
my $ct=1;
if ($opt_n){
	if ($opt_k){
		$ct=$opt_k+1;#we do this because we are  using 1 and not 0 as the first index into the array
	}else{ 
		my $info=`wc $opt_i -l`;chomp $info;
		$ct=int($info/$opt_B)+1;
	}
}else{
	if (!$opt_s){
#			system ("rm $opt_i\\_* -f\n");
		eval{
			my $targetfile=$opt_i;
			if ($opt_S){
				system ("sort $opt_i > tmp.sorted\n");
				$targetfile="tmp.sorted";
			}
			print "cmd:/usr/bin/split -l $opt_B -a 3 -d $targetfile $opt_i\\_ \n";
			system (" /usr/bin/split -l $opt_B -a 3 -d $targetfile $opt_i\\_ \n");  ### is this consistent across all machines???????
		};
		if ($?){
			die "Could not split $opt_i into blocks!\n";
		}
		$ct=`ls  | grep -P '$opt_i\_\\d\\d\\d\$' | wc -l`;chomp $ct;
	}elsif ($opt_k){
		$ct=$opt_k+1;
	}else{
		print "Running ls  | grep -P '$opt_i\_\\d\\d\\d\$' | wc -l\n";
		$ct=`ls  | grep -P '$opt_i\_\\d\\d\\d\$' | wc -l`;chomp $ct;
	}
	for (my $k=0;$k<$ct;$k++){
		$k=sprintf("%03d",$k);
		 my $ensembl=''; 
		 if ($opt_E){
			$ensembl=" --dbtype ensgene ";#die "using ensembl\n";
		}
        my $annot="$annovar --buildver $opt_o $opt_c --geneanno $ensembl --keepline --relpos --silent $opt_i\_$k $opt_d";
		if ($opt_v){
			$annot=~s/\-\-relpos/\-\-relpos \-\-novar /;
		}else{
#			 $annot.=";perl $bin/util/find_and_predict_NMD_indels.pl -i $opt_i\_$k.variant_function.collapsed ";
		}
		if (! -e "$opt_i\_$k.variant_function" && !$opt_n && !$opt_N){
			$id=printPBS("Annot\_$k",'8',$annot);
		}
	}
}
my $last_file="$opt_i\_".sprintf("%03d",$ct-1);#off by one bc we start at 0
my $lastfile_countline=`wc -l $last_file`;chomp $lastfile_countline;$lastfile_countline=~s/\s$last_file//g;
my $firstfile_headers_ct=`grep '^#' $opt_i\_000 | wc -l`;chomp $firstfile_headers_ct;$firstfile_headers_ct=~s/\s$firstfile_headers_ct//g;
print "There are $ct batches!!!\n";
print "The last file ($last_file) has $lastfile_countline and $firstfile_headers_ct\n";
#adjust the min memory and space requirements because you're now batching it
$filemem=-s "$opt_i\_000";$filemem=int($filemem/100000000);#in GB
#write on PBS job for each...
$minmem=2+$filemem;#this is the min amount of memory needed to run a small ANNOVAR database, larger files and databases will add more memory
if ($ct==1 && !$opt_S){
	unlink ("$opt_i\_000");
	system ("ln -s $opt_i $opt_i\_000\n");
}
my $fetch_script="$bin/../annofetch";print STDERR "Fetch script did not work..skipping ($fetch_script)" and $opt_z=0 if (! -e $fetch_script);
if ($opt_z && ! -e "$opt_i.flanks"){
	print STDERR "Running fetch script $fetch_script $opt_i $opt_i.flanks\n";
	$id=printPBS("fetch",'2',"$fetch_script $opt_i $opt_i.flanks");
	$opt_z="$opt_i.flanks";
}elsif ($opt_z && -e "$opt_i.flanks"){
	$opt_z="$opt_i.flanks";
	print "$opt_z already exists..using it\n";
}
my %db_size_hash;my $af;
if (!$opt_n){#run the databases in ANNOVAR
	for (my $k=0;$k<$ct;$k++){#make pbs bat files for all of the databases
		$k=sprintf("%03d",$k);
		if ($#files>-1 && !$opt_n){#runs if there are files and opt_n is not specified
			for (my $i=0;$i<=$#files;$i++){
				$files[$i]=~s/[\r\n]*$//;my ($db_dir,$db_file);
				print STDERR "working on $files[$i]\n";
				#find the size so that you can increase the memory on the grid
				if ($files[$i]=~/^\//){#if not absolute path, make it the absolute path
					$db_dir=`dirname $files[$i]`;chomp $db_dir;
					$db_file=`basename $files[$i]`;chomp $db_file;
				}else{
					$db_dir=$opt_d;
					$db_file=$files[$i];
					if (!-e "$opt_d/$files[$i]" && -e "./$files[$i]"){
						$db_dir="./";
					}
					if ($files[$i]=~/\_af$/){
						$af=$db_file;
						next;
					}
				}
				$flg{'sift'}++ if ($files[$i]=~/sift/i);
				$af++ and next if ($files[$i]=~/\_af$/);
				$convertible='hg19_hg18converted.txt' and next if ($files[$i]=~/hg19_hg18converted/);
				my ($int,$size);
				if (!exists $db_size_hash{"$db_dir/$db_file"}){
					$size= -s "$db_dir/$db_file";
					$int=int($size/100000000);#in GB
					$int+=$minmem; 
					$db_size_hash{"$db_dir/$db_file"}="$size,$int";
				}else{
					($size,$int)=split(",",$db_size_hash{"$db_dir/$db_file"});
				}
				if( ! -e "$db_dir/$db_file.idx" ){
					#try to index the database before beginning
					print "Trying to index file...$db_dir/$db_file\n";
					my $success=_index($db_dir,$db_file,$opt_m);
				}
				my $idxd=0;#indexed database exists
				if ( -e "$db_dir/$db_file.idx"){#indexing available so you do not need a lot of space
					$int="8";$idxd=1;
					$annovar .= " --batchsize 100k " if ($annovar!~/\-\-batchsize\s100k/i);
				}
				$db_file=~s/$opt_o\_//g;$db_file=~s/.txt//g;#$files[$i]=`basename $files[$i]`;chomp $files[$i];
				my $addon;
				if ($opt_v){$addon=' --novar ';}#if your input file has no variants and you just want to annotate regions based on coordinates
				my $toDo="$annovar $opt_c --buildver $opt_o --dbtype $db_file $addon --regionanno --keepline --silent $opt_i\_$k $db_dir";
				
				print "PBS$k\_$i...($size)$int\n";
				if (($db_file=~/(1000g|snp1[23]\d|avsift|ljb|vcf|cg69|siftv63|hapmap)/i)){#large files with index
					$toDo="$annovar --buildver $opt_o --dbtype $db_file --filter --keepline $addon --silent $opt_i\_$k $db_dir";
					$int=8;#we do this because it is filter option and memory should not be an issue
				}elsif ($db_file=~/(cosmicdb|SomamiR_targets|SomamiR)/ && $use_strict){
					$toDo="$annovar --buildver $opt_o --dbtype generic --genericdbfile $opt_o\_$db_file.txt --filter --keepline $addon --silent $opt_i\_$k $db_dir";
				}elsif ($db_file=~/(hg19_.*hom_only|hg19_CG-.*AllFreq)/ &&  $use_strict){
					$toDo="$annovar --buildver $opt_o --dbtype generic --genericdbfile $db_file --filter --keepline $addon --silent $opt_i\_$k $db_dir";
					$int=8
				}elsif ($size > 1_000_000 && !$idxd){
				}
				my $id=printPBS ("$k\_$i",$int,$toDo);
			}#ends foreach file
		}#end if statement else should be $#files && $opt_n}
#		}elsif ($#files && $opt_n){
#			print "[INFO] Skipping running databases for ".join (",",@files)."\n";
#		}#ends if/else
	}#end foreach k loop that sets up all the interim files with the database files;

	$notDone=1;
	$notDone=monitor ($jobs_qry) unless ($jobs_qry !~/\d+/);
	if ($opt_g){$notDone=0;}#not monitoring anymore
	if ($notDone){
		die "Your grid jobs exceeded the amount of time or did not get submitted correctly to the grid...monitor your jobs manually and concatenate your results\n" ;
	}
}
#check to see which ones ran properly
# Max processes for parallel download ($MAX_PROCESSES)
my %err;#logs errors
my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
print "About to fork processes\n";
my %donotaggregate;
CHUNK: foreach (my $k=0;$k<$ct;$k++){
	$k=sprintf("%03d",$k);
	my @outputfiles;
	$pm->start and next; # do the fork
	my $annovar_output="./$opt_i\_$k.variant_function";
	for (my $i=0;$i<=$#files;$i++){
		my $file;
		#locally run
		print "cannot be null...at line ".__LINE__ ."\n" and next if ($files[$i] eq '' || $files[$i]=~/^\s+$/);
		chomp $files[$i];
		my ($dirname,$filename)=(`dirname $files[$i]`,`basename $files[$i]`);chomp $dirname;chomp $filename;
		$filename=~s/\r//g;
		until ($filename!~/(.*)\.txt/){
			$filename=$1;
		}
		$flg{'sift'}++ if ($files[$i]=~/sift/);
		$af++ and next if ($files[$i]=~/hg19_af$/);#renamed the af
		next if ($files[$i]=~/$convertible$/);
		print "Running system: ls | grep $opt_i\_$k\.$filename\$ \n";
		$file=`ls | grep $opt_i\_$k\.$filename\$`;chomp $file; 
		print "looking for $filename...found($file)\n";
		
		if (!$file|| -z $file){
			print "Running system(2nd attempt): ls | grep  $opt_i\_$k\.$filename | grep -e filtered -e exact | grep -ve tmp\n";
			$file=`ls | grep  $opt_i\_$k\.$filename | grep filtered -e exact |grep -ve tmp`;chomp $file;
			print "\tworking for $file\n";
		}
		if (!$file){	
			$file="$opt_i\_$k.$filename";
			print STDERR "======================\nCould not find a file for Database: $filename and Outputfile: $opt_i\_$k.$filename...\n======================\n"; 
			my $tmpname=$filename;
			if ($filename=~/(hg19\_|mm\d\_)(\w+)/){
				$tmpname=$2;
			}else{
				die "$filename does notmatch regex\n";
			}
			my $success=_makeERRFile($file,$tmpname);
			if ($success){
				$err{$annovar_output}.="Internal Error: Could not find the db run for CHUNK $k and db $file\n";
			}
			#next CHUNK;
		}#CHUNKYMONKEY
		
		if ($file=~/filtered$/ && !-z $file){
			#sort it and then output only the score
			$file=_sortFilterFile($file);
		}
		push (@outputfiles,$file);
	}
	print STDERR  "[INFO] Finished collecting data on databases across ".sprintf ("%d",($k-1))."...:".join (",",@outputfiles)."\n";
	
	if (-e "$annovar_output.collapsed" || $opt_N){
	}elsif (-e $annovar_output ){
		if ($opt_l && -e "$annovar_output.collapsed"){#collapse the exonic and variant file beginning at $opt_l
			collapseANNOVAR($annovar_output,$opt_l);
		}elsif ($opt_n && -e "$annovar_output.collapsed"){ #did not run annot and the file already exists
			#do nothing, hopefully this is your completed run
		}elsif ($opt_v && !-e "$annovar_output.collapsed"){
			open (COLLAPSED, ">$annovar_output.collapsed" ) or die "Cannot open $annovar_output.collapsed\n";
			my $annot_ver=($opt_E)?"Ensembl v.63":"NCBI v.37";
			print COLLAPSED join ("\t","#ANNOVAR annot $annot_ver","Gene","Rel pos to feature(s)","Gene Ori","Chr","Query Start","Query End","Allele1","Allele2","Comment")."\n" ;
			close COLLAPSED;
			system ( "cat $annovar_output >>$annovar_output.collapsed\n");
		}else{#ran annot and collapse
			collapseANNOVAR($annovar_output,0);
			system ("rm $annovar_output\.*_variant_function -f \n") if (!$opt_r);
		}
	}else{
		print STDERR "$annovar_output does not exist and -N and $annovar_output.collapsed does not exist\n";
		my $success=_makeERRFile("$annovar_output.collapsed",join ("\t",'ANNOVAR annot','Gene','Rel pos to feature(s)','Gene Ori','Chr','Query Start','Query End','Allele1','Allele2','Comment'));
		if ($success){
			$err{$annovar_output}.="Internal Error: Could not find the db run for CHUNK $k and gene annot file: $annovar_output.collapsed\n";
		}
	}
	
	print STDERR "[INFO] About to paste all samples for $k\n";
	if ($#outputfiles>-1){
		eval{
			my $cmd="paste ".join (" ",@outputfiles) ;
			if (!$opt_N){
				$cmd.=" $annovar_output.collapsed > $out\_$k.output";
			}else{
				$cmd.=" > $out\_$k.output";
			}
			print "$cmd\n";
			system ("$cmd\n");
		};
		if ($@){
			die "could not join files together\n";
		}elsif (!$opt_r){
			eval{
				system  ("rm $annovar_output.collapsed ". join (" ",@outputfiles)."\n");
				system  ("rm PBS* -f ");
			};
			if ($@){
				warn "could not rm intermediate files\n";
			}
		}
	}else{
		eval{
			system ("cp $annovar_output.collapsed $out\_$k.output\n");
		};
		if ($@){
			die "could not cp $annovar_output.collapsed\n";
		}elsif (!$opt_r){
			system ("rm $annovar_output.collapsed\n");
		}
	}
	if ($opt_I || (-e "config.ini" && `grep 'runIndelNMD=on' config.ini | wc -l`>=1)){
		system ("perl $bin/../find_and_predict_NMD_indels.pl -i $out\_$k.output -o $out\_$k.PredIndel \n");
		print "Running $bin/../find_and_predict_NMD_indels.pl -i $out\_$k.output -o $out\_$k.PredIndel\n";
		if (-e "$out\_$k.PredIndel" && !-z "$out\_$k.PredIndel"){
			if (!$opt_r){eval{system ("mv $out\_$k.PredIndel $out\_$k.output\n");};
				if ($?){
					warn "could not move $out\_$k.PredIndel to $out\_$k.output(ok)\n"; 
				}
			}else{
				system("cp $out\_$k.PredIndel $out\_$k.output\n");
			}
		}
	}
	print STDERR "[INFO]  Done Pasting for $k\n";
	$pm->finish; # do the exit in the child process
}
print STDERR "[INFO] Waiting for children processes to complete\n";
$pm->wait_all_children;
if (keys %err){
	print "The following errors occurred during aggregation:\n";
	foreach my $errmsg (keys %err){
		print "$errmsg\n";
	}
}
my @append=split("\n",`ls $out\_* | grep output\$`);
print STDERR "($out.output)working on appending the following:\n\t".join("\n\t",@append)."\n";
if ((-e "$out.output" && -z "$out.output") || !-e "$out.output"){
	print STDERR "Working to get headers\n";
	system ("grep '^#' -m 1 $append[0]>$out.output\n");
	if ($opt_R){
		system ("perl $bin/rename_headers.pl $out.output\n");
	}
	system ("cp $out.output $out.headers.output\n");
}
#die "Concatenating...". join (",",@append)."\n";;
my $err_flag;
for (my $d=0;$d<=$#append;$d++){
	chomp $append[$d];
	eval {
		system ("grep -ve '^#' $append[$d] >>$out.output\n");
		#time trial to see which method is quicker
		#system ("perl -pi -e '$_ = "" if ($. == 1);' $append[$d]\n");
		#system ("cat $append[$d] >> $out\n");
	};
	if ($?){
		warn "could not append $append[$d] to $out.output\n\tRan:grep -ve '^#' $append[$d] >>$out.output\n$?";
		$err_flag.="\n\t[INTERNAL ERR] The number of input files does not equal output files. Lines ". (($d*$opt_B)+1) . "-". (($d+1)*$opt_B)." ($d and $opt_B) myfailed.  Please resubmit entire file or selected lines to AVIA.\n";
	}else{
		if (!$opt_r){
			eval{
				system ("rm $append[$d] -f \n") ;
			};
			if ($?){
				warn "could not remove intermediate file $append[$d] at LN".__LINE__."\n";
			}
		}
	}
}
if ($opt_N){
	eval{
		system ("echo '#Comments' > test.collapse;cat $opt_i>>test.collapse;paste $out.output test.collapse>> test.output;mv test.output $out.output\n");
	};
	if ($?){
		warn "Could not append the input file to the outputfile\n";
	}
}
cleanup($ct);

print STDERR "Done!\n";
$af=" ANNOVAR.input.hg19_af " if ($af);#assume that this is from the AVIA pipeline and no files have been renamed
$af='' if (!-e $af || -z $af);
my $addon='';
if ($opt_z && !$err_flag){
	open (FLANK,"<$opt_z") or die "cannot open $opt_z\n";
	my @id=<FLANK>;
	if ($id[0]=~/ERR/ && ( -e "cvrt2anvr.stderr.log" || `head $opt_i -n1`=~/^#/)){ 
		$id[0]=~s/ERR\s{0,}\n/#FlankingSequence\n/;
	}else{
		unshift(@id,"#FlankingSequence\n");
	}
	close FLANK;
	open (FLANKWRITE,">$opt_z.wHeaders") or die "cannot open $opt_z.wHeaders";
	print FLANKWRITE join ("",@id);
	close FLANKWRITE; 
	$addon="$opt_z.wHeaders ";
#	system ("paste $opt_z.wHeaders $af $out.output> $out.output1;mv $out.output1 $out.output;\n");
}
if ($af){
	$addon.="$af ";
	
}
if (-e "config.ini" && -e "$convertible" && !-z "$convertible"){
	$addon.="$convertible ";
}
if ($addon){
	system("paste $addon $out.output> $out.output1;mv $out.output1 $out.output;\n");
}

if ($err_flag){	
	open (FINAL,">>$out.output") or die "Cannot open $out.output for appending error message at LN".__LINE__." in $0\n";
	print FINAL "\n\nThe following errors occurred during processing:\n$err_flag";
	close FINAL;
}
if ($opt_R || (-e "config.ini" && `grep '^name_column_header=1' config.ini | wc -l`==1)){
	my ($input_header_name,$id);
	if ($id=`grep ^file= config.ini`){
		$id=`basename $id`;
		$id=~s/(^\d+\.|\.txt|\.vcf)//g;
	}elsif ($id=`grep '^input_fullpath=' config.ini`){
		$id=`basename $id`;
		$id=~s/(\.txt|\.vcf)//g;
	}else{                                                                                
		$id=$opt_i;                                                                   
		$id=~s/(hg19\_|\_Blood|\_tumor)//i;                                           
		$id=~s/(^\d+\.|\.txt|\.vcf)//g;                                               
	}       
	$id=~s/^(\d{14}\.)//g;
	if ($id){
		chomp $id;
		eval{
#			system (" awk -F '^' '{print \"$id\",\$1}' OFS='\\t' $out.output > .$out.output\n");
			system ("awk -F '^' '{if ( NR==1 )  print \"#SampleName\",\$1  ; else print \"$id\",\$1  }' OFS='\\t' $out.output > .$out.output\n");
		};
		if ($@){
			print STDERR "ERROR Could not append to the output your input filename\n";
		}
		if ($opt_R && !-z ".$out.output"){
			 system ("mv .$out.output $out.output\n");
		}
	}
}

system ("touch DONE.$opt_i.anvr.wrpr\n") if ((-s "$out.output") > 101);
exit();
sub _makeERRFile{
	my $xfn=shift;
	my $xhead=shift;
	my $xfct;
	if ($xfn=~/\_(\d{3})[\_\.]/){
		$xfct=sprintf("%d",$1);
	}
	my $flag;
	if ($xfct==$ct-1){
		$flag=$lastfile_countline-1;
	}else{
		$flag=$opt_B-1;#print $flag;
	}
	if ($xfct==0){#subtract the headers starting with '#' too!
		$flag=$flag-$firstfile_headers_ct;
	}
	print STDERR "Could not run  _makeERRFile for an empty filename\n" and return 0 if (!$xfn);
	print STDERR "Running _makeERRFile for $xfn($flag)\n";
	open (ERRDFN,">$xfn") or die "Cannot open file for err\n";
	print ERRDFN "#$xhead\n";
	for (my $i=0;$i<=$flag;$i++){
		print ERRDFN "INTERNAL ERR \n";
	}
	close ERRDFN;
	return 1;
}
sub cleanup{
	my $ct=shift;
	return if (!$opt_r);
	#this subroutine moves the PBS jobs and intermediate files into a subdirectory
	my $dir=`date "+%Y%m%d"`;chomp $dir;
	system ("mkdir $dir\n");
	system ("mv PBS*bat PBS*bat.o* *log *valid_input *notexon.input $dir/ \n"); 
#	if ($ct>0){
#		for (my $i=0;$i<=$ct;$i++){
#			$i=sprintf("%03d",$i);
#			system ("rm $opt_i\_$i $dir\n");
#		}
#	}
	return;
}
sub monitor{
	my $targets=shift;
	if ($opt_g){
		return 0;
	}
	my $timeout=`date "+%d%H%M%S"`;chomp $timeout;$timeout+=$opt_T;#1hour
	my $done=0;
	until (`$targets` eq '' || $done){
		sleep (60);
		my $curr_time=`date "+%d%H%M%S"`;chomp $curr_time;
		if ( $curr_time>$timeout){$done=1;}
	
	};
	return $done;
}
sub _populate_rank{
	#rank the impacts
	%rank=('3PRIME_UTR'=>3,
'5PRIME_UTR'=>3,
'NMD_TRANSCRIPT'=>1,
'SPLICE_SITE'=>2,
'CODING_UNKNOWN'=>1,
'DOWNSTREAM'=>5,
'ESSENTIAL_SPLICE_SITE'=>1,
'FRAMESHIFT_CODING'=>1,
'INTERGENIC'=>5,
'INTRONIC'=>4,
'NON_SYNONYMOUS_CODING'=>1,
'SYNONYMOUS_CODING'=>2,
'REGULATORY_REGION'=>3,
'WITHIN_NON_CODING_GENE'=>2,
'STOP_GAINED'=>1,
'UPSTREAM'=>4
	);
}
sub printPBS{
	my $name=shift; my $memsize=shift; my $cmd=shift;
	open (OUT ,">$cwd/PBS.$opt_i.$name.bat") or die "Cannot open $cwd/PBS.$name.bat for writing\n";
	print OUT  "#!/bin/sh\n";
	print OUT  "#PBS -j oe -l arch=linux-x86_64,pvmem=$memsize\GB,mem=$memsize\GB\n";
	print OUT  "#PBS -l cput=10:000:00,pcput=10:00:00,walltime=10:00:00\n";#10 hours
	print OUT  "cd $cwd\n";
	print OUT  "$cmd \n";
	close OUT ;
	my $xid=' ';
	#launch the job and return the id
	if (!$opt_g){
		my $jobs=`qsub PBS.$opt_i.$name.bat`;
		chomp $jobs;
		if ($jobs=~/(\d+).abcc1.ncifcrf.gov/i){
			$xid=$1;$jobs_qry.=" -e $xid ";
		}else{
			print STDERR "[ERR] Could not submit PBS.$name.bat to the grid\n";
		}
		push(@gids,$xid);
	}else{
		eval{
			system ("chmod +x PBS.$opt_i.$name.bat;./PBS.$opt_i.$name.bat\n");
		};
		if ($?){
			warn "Error processing chmod and running command on PBS.$opt_i.$name.bat in ". `pwd`;
		}
		push (@gids,'1');
		$xid='1';
	}
	return $xid;
}
sub checkPBSAccess{
	if ($opt_g){#force local run
		$opt_g=1;return;
	}
	print STDERR "checking PBS...\n";
	eval {
		system ("qstat\n");#checking if qstat errors out, if so , cannot run on grid
	};
	if ($?){
		$opt_g=1;
		print STDERR "[WARN] You cannot use the grid to run your jobs (system err)$?\n";
	}else{
		print STDERR "[INFO] Submission to grid is allowed\n";
	}
	return;
}
