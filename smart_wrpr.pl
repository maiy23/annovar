#!/usr/bin/perl
use strict;
use FindBin;
use Data::Dumper;
use Getopt::Std;
use lib '/users/abcc/vuonghm/scripts.dir/perl_modules/';
umask(0000);
my $promptline=join (" ",@ARGV);
use vars qw($opt_o $opt_f $opt_a $opt_t $opt_g $opt_A $opt_G $opt_v $opt_w $opt_O $opt_V $opt_F $opt_N $opt_j);
getopts("a:t:o:f:A:G:O:j:gvwVFN");#database and avia id, respectively
=head
This is the wrapper for several wrappers which is specified by the tool $opt_t.  This runs funseq, and outputs into FDI output
=cut
my $bin='/bioinfoC/hue/scripts.dir/';chomp $bin;
my $usage=qq(
	$0 
	[REQUIRED]
	-a <fullpath to folder name>
		#it is assumed that the directory will have a file called <opt_a> (default from avia pipeline)
		#if not then specify the input filename with -f option
	-t <tool to run>
		Accepted names:
		--funseq
		--intersectBed
		--merge
		--uniprot
	[OPTIONAL]
	-o <options list for tool>
	-f <input filename>
	-g <if specified; run on grid>
	-A <filename of list of databases>
	-G <groupBy options>  #only used if intersectBed is run;otherwise ignored
		#DEFAULT settings:  	-grp 1,2,3,4,5 -opCol 8 -op collapse
		#this option collapses by the first -grp columns and collapses based on -opCol
		#see bedtools for more options
	-v <if specified, then verbose>
	-w <if specified, run from web>
	-O <output filename base>
	-V/F <in AVIA/FDI format>
		--the '-F' option works with funseq only
	-N <no headers>
	-j <for FunSeq2 only, specifies the name of the output directory>

);
my $fs2outputdir=(defined $opt_j)?$opt_j:'out';
die $usage if (!defined $opt_a || !defined $opt_t);
my $id=`basename $opt_a`;chomp $id;
my $dirname;
if ($opt_a!~/^\//){
	$dirname='/bioinfoC/AVA/FDI/'.$id;
}elsif (!defined $opt_f){
	$dirname=$opt_a;
}else{
	$dirname=$opt_a;
}

$opt_f||="$id";
if ($opt_w){$opt_f="ANNOVAR.input";$opt_O="$id";}#run for all
if (!$opt_v){
	open (LOG ,">>$dirname/smrt_wrpr.log")  or die "Cannot open log file\n";
	open (STDERR, ">&LOG");
}
print STDERR "[INFO] Starting $0 at ". `date`;
print STDERR "[COMMAND LINE] $0 $promptline\n";
################ Set up the tool commands here ##################3
my $funseq_exe="/bioinfoC/hue/scripts.dir/perl_modules/funseq-0.1/scripts/funseq.sh";# -f All.input.txt -f 2 -inf bed -outf bed -nc
my $funseq2_exe="/SeqIdx/FunSeq2/funseq2-1.0/run.sh ";# -f All.input.txt -f 2 -inf bed -outf bed -nc
my $intersectBed_exe='/opt/nasapps/development/bedtools/2.17.0/bin/intersectBed';
my $wrapper='/bioinfoC/hue/scripts.dir/parseData/mergeBedByNcols.pl -u ';
my $groupBy_exe='/opt/nasapps/development/bedtools/2.17.0/bin/groupBy';
##################################################################
my $DB_LOC='/SeqIdx/annovardb/humandb';
chdir($dirname);
print STDERR "I'm currently in $dirname and using $opt_f\n" ;my $out;
if ($opt_t=~/funseq2/i){
	#run funseq
	$opt_o||=' -m 2 -inf bed -outf bed -db';
	$out=(defined $opt_O)?$opt_O:"$opt_f.hg19_FunSeq2";
	print "FunSeq2 out file: $out\n";
	if (! -e "$out"){
		my $cmd;$id=($id=~/^([\w\-]+)\./)?$1:$id;
		if (!-e "$fs2outputdir/Output.bed" && !-e "$fs2outputdir/Output.indel.bed"){$cmd="$funseq2_exe $opt_o -f $opt_f -o $fs2outputdir\n";monitor_grid($cmd) ;}
		if (-e "$fs2outputdir/Output.bed"  && -e "$fs2outputdir/Output.indel.bed"){
			mergeFS2 ("$fs2outputdir/Output.bed","$fs2outputdir/Output.indel.bed");
		}elsif (-e "$fs2outputdir/Output.bed" ){
			mergeFS2 ("$fs2outputdir/Output.bed");
		}elsif (-e "$fs2outputdir/Output.indel.bed" ){
			mergeFS2 ("$fs2outputdir/Output.indel.bed");
		}else{
			die "Your file Output.bed did not generate\n";
		}
	}else{
		print STDERR "$out already exists...just printing:\n";
		printout("$out");
	}
}elsif ($opt_t=~/funseq/i){
	#run funseq
	$opt_o||=' -m 2 -inf bed -outf bed ';
	$out=(defined $opt_O)?$opt_O:"$opt_f.hg19_FunSeq";
	if (! -e "$out"){
		my $cmd;$id=($id=~/^([\w\-]+)\./)?$1:$id;
		if (!-e "out/$id.FUNSEQ.bed"){$cmd="$funseq_exe $opt_o -f $opt_f \n";monitor_grid($cmd) ;}
		my $tmp=`echo $opt_f  |cut -f1 -d '.'`;chomp $tmp;
		if (-e "out/$tmp.FunSEQ.bed"){
			merge ("out/$tmp.FunSEQ.bed","$opt_f","$out");
		}else{
			die "Your file out/$opt_f.FunSEQ.bed did not generate\n";
		}
	}else{
		print STDERR "$out already exists...just printing:\n";
		printout("$out");
	}
}elsif ($opt_t=~/intersectBed/i){
	#first things first make sure that your input file is not junky and begins with 'chr'
	#for the encode data especially, it will not intersect if it doesn't have the 'chr' in the front!
	die "Your specified file $opt_f does not exist\n" if (! -e $opt_f);
	my $newname="$opt_f.xbed";
	my @dbnames;
	my $input_has_chr=`grep  -P '^chr[\\dXYM]' $opt_f | wc -l`;chomp $input_has_chr;
	#on the command line, you should specify which file to use (use the -o option)
	# /opt/nasapps/development/bedtools/2.17.0/bin/intersectBed -a foo -b /SeqIdx/annovardb/humandb/hg19_wgEncodeCrgMapabilityAlign100mer.txt -wb -wa > moocowplease"
	if (!defined $opt_o && !defined $opt_A){
		die "You must specify -b option for intersectBed (this is the db file you wish to intersect your file with!  Use -A <fn> or -o <list of options including -b nameOfDB>)\n";
	}else{
		my $addon_opts=" -wb -wa ";
		if ($opt_o=~/\-b\s+(\S+)/){
			my $fullpath=$1;
			my $dbname;
			if (! -e "$fullpath"){#not in the current working directory so add it
				$dbname=$fullpath;
				$fullpath="$DB_LOC/$fullpath";
				$opt_o=~s/\-b\s+\S+//;
			}else{
				$dbname=`basename $fullpath`;chomp $dbname;
			}
			push (@dbnames,$dbname);
		}elsif ($opt_A){
			open (DBS,"<$opt_A") || die "Cannot open file\n";
			@dbnames=<DBS>;
			close DBS;
			print STDERR "DATABASES:\n".join ("",@dbnames) ;
		}else{
			die "You must specify an input filename in the -o option of this script to be passed to the intersectBed (-b DBFILENAME)\n";
		}
		if ($opt_o!~/-wb/){
			$opt_o.=" -wb ";
		}
		if ($opt_o!~/-wa/){
			$opt_o.=" -wa ";
		}
		foreach my $dbname (@dbnames){
			chomp $dbname;
			$out=(defined $opt_O)?$opt_O:$opt_f;
			my $outname=$out.".".`basename $dbname | sed 's/\.txt\$//g' `;chomp $outname;
			# print "writing to $out($outname)\n";<STDIN>;
			if (!-e $dbname && ! -e "$DB_LOC/$dbname"){
				print STDERR "Could not find your database $dbname\n" ;next;
			}elsif (!-e $dbname){
				$dbname="$DB_LOC/$dbname";
			}
			my ($success,$features,$startsWith)=checkDB($dbname);
			next if (!$success && ! $features);
			print STDERR "[INFO] The column that we want to output is $success and the features are $features;$startsWith\n" ;
			my $xfile;
			if ($startsWith){#if db with chr 
				$xfile="$out\_wchr.tmp";
				if (! -e $xfile){
					if ($input_has_chr){#&& input starts with chr
						print STDERR "TEST1($input_has_chr)\n" ;
						_system  ("awk '/^(chr[1-9XYMT]*)/' $opt_f  | awk '!/^(chr|#)/{sub(/^/, \"chr\")}; 1'|cut -f1-5 > $xfile\n") ;
					}else{
						print STDERR "TEST2($input_has_chr)\n" ;
						_system ("awk '!/^(chr|#)/{sub(/^/, \"chr\")}; 1' $opt_f |cut -f1-5 > $xfile\n") ;
					}
				}
			}else{#database does not start with chr
				$xfile="$out\_nochr.tmp";
				if (! -e $xfile){
					if ($input_has_chr){#input starts with chr
						print STDERR "TEST3 ($input_has_chr)\n" ;
						_system ("awk '/^(chr[1-9XYMT]*)/' $opt_f  | awk '!/^(#)/{sub(/^chr/, \"\")}; 1'|cut -f1-5 > $xfile\n") ;
					}else{
						print STDERR "TEST4 ($input_has_chr)\n" ;
						_system ("awk '!/^(#)/{sub(/^chr/, \"\")}; 1' $opt_f|cut -f1-5 > $xfile\n");
					}
				}
			}	
			 monitor_grid ("$intersectBed_exe $opt_o -b $dbname -a $xfile> $outname\n");
			# print "$intersectBed_exe $opt_o -b $dbname -a $xfile> $outname\n";<STDIN>;
			if (!$opt_G){
				$opt_G=" -g $features -c ". ($success+5)." -o collapse";#default if input is 5 columns and the 4th column in db is score
			}
			if (-e "$outname" && ! -z "$outname"){
				# print "$outname is not equal to zero \n". `ls $outname -al `;exit;
				monitor_grid  ("$groupBy_exe $opt_G -i $outname > $outname.collapsed\n");
				if ( -e "$outname.collapsed" && ! -z "$outname.collapsed"){
					my $out=($opt_w)?' -F ':'';
					merge ("$outname.collapsed","$opt_f","$outname",substr($features,-1));
				}
			}elsif (-z $outname){print "This is the end\n";
				_system ("mv $outname $outname.tmp\n");
				merge ("$outname.tmp","$opt_f","$outname",substr($features,-1));
			}
			# system ("rm *tmp -f ;mv $outname.collapsed $outname\n");
			print STDERR "Completed with $dbname\n" ;
		}

	}
	system ("chmod 777 $opt_a -R 2>/dev/null\n");
}elsif ($opt_t=~/uniprot/){
	$opt_O||="$opt_f.uniprot_results.txt";
	_system( "perl $bin/parseData/convert4Uniprot.pl -a $opt_f -m \n");#collapsed annovar file
	_system( "php /mnt/webrepo/fr-s-abcc-bda-0/htdocs/mongoDb/mongo_find.php  -t tab -v  $opt_f.prot |grep -ve '{' -ve '^\$' | sed 's/^ //'> $opt_f.uniprot.txt\n");
	merge("$opt_f.uniprot.txt","$opt_f.prot.mapping","$opt_O"," 1 -c 7 -k ")
}else{

}
sub _system {
	my $cmd=shift;
	print STDERR "[EXE] $cmd" ;
	eval{
		system ("$cmd");
	};
	if ($? || $@){
		warn "[ERR] Could not run $cmd\n\n" ;
	}else{
		print "[INFO] Successful cmd run ($cmd)\n\n" ;
	}
}
sub merge{
	my ($file1,$file2,$outputname,$field)=(@_);
	my $addon='';
	if ($field){
		$field=" -n $field ";
	}else{
		$field=' ';
	}
	
	if ($opt_V){
		$field.=" -A ";
	}
	if ($opt_N){
		$field.=" -N ";
	}
	if ($opt_F){
		$field.=" -F ";
	}
	print "About to run perl $wrapper -f $file1 -s $file2 -o $outputname $addon $field\n";
	_system ("perl $wrapper -f $file1 -s $file2 -o $outputname $addon $field\n");
	return;

}
sub printout{
	my $xfile=shift;
	open (XFILE,"<$xfile") or die "Cannot open $xfile\n";
	while (<XFILE>){
		print $_;
	}
	close XFILE;
}
sub monitor_grid{
	my $cmd=shift;
	open (BAT,">$$.bat" ) or die "Cannot open $$.bat\n";
	print BAT "#PBS /bin/bash\nPBS_O_WORKDIR=$dirname\n#PBS -j oe -e $opt_a/$$\_stderrout\ncd \${PBS_O_WORKDIR}\numask 000\n";
	print BAT "$cmd";
	close BAT;
	system ("chmod 777 $$.bat 2>/dev/null\n");
	print  "Running command $cmd in runAndMonitor\n";
	if ($opt_g){
		my $grid_id=`qsub $$.bat`;chomp $grid_id;
		print  "In runAndMonitor...$cmd($$.bat)...grid_id=$grid_id\n" ;
		if ($grid_id=~/(\d{6,})/){
			my $targets="qstat | grep $1 ";
			my $timeout=3600;#1hour=6000
			my $done=1;
			until (`$targets` eq '' || !$done){
				sleep (15);
				my $elapsed= (time() - $^T);
				if ( $elapsed > $timeout){$done=0;print LOG "Timed out\n";}
			
			};

			return $done ;
		}else{
			print STDERR "Could not submit to grid ($$.bat) for  $cmd and ($grid_id)\n" ;
			die "Could not submit to grid ($$.bat) for  $cmd and ($grid_id)\n";
		}
	}else{
		print STDERR "Running locally $$.bat\n" ;
		system ("./$$.bat\n");
	}
	return 0;
}
sub checkDB{
	my $dbname=shift;
	if (-e $dbname){
		my $chr=`head -n1 $dbname`;
		my $file=`basename $dbname`;chomp $file;$file=~s/(hg19\_|.txt)//g;
		#first find out which column is the output column for -c option
		print STDERR "looking for $file in $DB_LOC/lkup.txt\n" ;
		my @idx=split("\t",`grep $file $DB_LOC/lkup.txt`);
		$idx[4]++ if ($idx[4]);##lkup.txt if offset by one (starts with 0)
		#column 4 is the output column
		my $startsCHR=0;
		if ($chr=~/^1/){
			#find out if there are alleles
		}elsif ($chr=~/^chr/){
			#assume chr
			$startsCHR=1;
		}elsif ($chr=~/^\d+/){
			#cannot use
			return 0;
		}else{
			#don't know what else to expect so we ignore
			return 0;
		}
		my @row_elems=split("\t",$chr);
		if ($row_elems[3]=~/^[ATCG]*$/i){
			$idx[4]||=6;
			return ($idx[4],'1,2,3,4,5',$startsCHR);#yes it has alleles
		}else{
			$idx[4]||=4;
			return ($idx[4],'1,2,3',$startsCHR);#no it does not have alleles
		}

	}else{
		print STDERR "[WARN] $dbname did not exist...Skipping\n" ;
		return 0;
	}
}

sub mergeFS2{
	my (%hash,%headers);
	print "Files:",join ("," ,@_);print "\n";
	foreach my $file (@_){
		print "working on $file\n";
		print "File $file doesn't exist!??" and next if (! -e $file);
		open (FH,"<$file") or die "Cannot open $file\n";
		my $headerFound=0;
		while (<FH>){
			my @arr=split("\t",$_);
			chomp $arr[$#arr];
			my @columns=split(";",$arr[6]);
			if ($headerFound){
				my $idx=join("\t",@arr[0..4]);#chr start stop ref var
				$idx=~s/^chr//g;
				if (exists $hash{$idx}){
					print "$idx already exists" ;
					# die "$idx already exists ". Dumper (\%{$hash{$idx}}) ."\n". join ("\n",@columns);
				}
				foreach (my $i=0;$i<=$#columns;$i++){
					$hash{$idx}{$headers{$i}}=$columns[$i];
				}
			}else{
				if ($_=~/^#/){
					print "Header found\n";
					$headerFound++;
					foreach (my $i=0;$i<=$#columns;$i++){
						$headers{$i}=$columns[$i];
						$headers{"$columns[$i]"}=1;
					}
				}else{
					die "Could not find header on first line\n";
				}
			}
		}
		close FH;
		$headerFound=0;
	}
	# print "[DEBUG] waiting for files to be done reading...\n" and <STDIN>;
	open (FILE,"<$opt_f") or die "Opening $opt_f...\n";
	print "opening $opt_f\n";
	open (OUT,">$out") or die "Cannot write to $out\n";
	print "writing to $out\n";
	
	open (HEADER,">FunSeq2.headers.out") or die "Cannot open FunSeq2.headers.out\n";
	my $feat;
	foreach my $id (sort keys %headers){
		next if ($id=~/^\d+$/);
		$feat.= "$id;";
	}
	chop $feat;
	print HEADER "$feat\n";
	close HEADER;
	if (!$opt_F){
		print OUT "#";
		foreach my $id (sort keys %headers){
			next if ($id=~/^\d+$/);
			print OUT "$id;";
		}
		print OUT "\n";
	}
	my $printer='';
	while (<FILE>){
		next if ($_=~/^#/);
		my @arr=split("\t",$_);chomp $arr[$#arr];
		my $idx=join("\t",@arr[0..4]);
		$idx=~s/^chr//;
		if ($opt_F){
			my $j=$idx;
			$j=~s/\t/:/g;
			$printer="$j\t";
		}else{
			$printer='';
		}
		foreach my $id (sort keys %headers){
			next if ($id=~/^\d+$/);
			if (exists $hash{$idx}{$id}){
				$printer.= "$hash{$idx}{$id};";
			}else{
				$printer.=".;";
			}
		}
		$printer=~s/;$//;
		print OUT "$printer\n";
	}

}