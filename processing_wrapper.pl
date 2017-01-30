#!/usr/local/bin/perl
=head v2.0
This version is uses the fact the website accepts gzip, bzip, tar, tgz file formats for upload
So we must sanitize the inputs before moving to our work directory  (for ia_upload_only)
No longer accepts lmt inputs for checks!

v3.0
addition of viz options in the mix
=cut
use strict;
use FindBin;
use lib "$FindBin::Bin/";
use lib "$FindBin::Bin/perl_modules/Config-IniFiles-2.38";
use lib "$FindBin::Bin/perl_modules";
use lib "/bioinfoC/AVA/prod_scripts/util/perl_modules";
use Mail::Sendmail;
use Config::IniFiles;
use Getopt::Std;
use parse;
use Data::Dumper;
umask(0000);
use vars qw( $opt_o $opt_f $opt_c $opt_t $opt_w );
getopts("f:c:t:o:w");
########### FIND EXECUTABLES FOR RUNNING PIPELINE (make sure they exist on WEB SERVER!!!) ############
my $tar=`which tar`;chomp $tar;
my $gunzip=`which gunzip`;chomp $gunzip;
my $bzip=`which bunzip2`;chomp $bzip;
my $unzip=`which unzip`;chomp $unzip;
my $ADMIN="vuonghm\@mail.nih.gov";
my $perl=`which perl`;chomp $perl  ; 
my $qsub=`which qsub`;chomp $qsub;
my $_rename=0;	my @files;
my $refver;#this can be changed during the course of this script where as $cfg{'USER'}{'refver'} is NOT, preserves original request info
my $bin=`dirname $0`;chomp $bin;
my $dbinfofile='/SeqIdx/annovardb/XXXX/database.info';
umask(0000);

if (!$qsub && ! -e "/usr/local/torque/bin/qsub") {
	 die "could not run this executable at LN $0 ".__LINE__." ($qsub) and /usr/local/torque/bin/qsub does not exist on server \n";
}elsif (!$qsub){
	$qsub="/usr/bin/qsub";
}
my $usage=qq(
	$0 -f <web_input> -o <config_output>
);
#SAFE_PROCESS ("source /users/abcc/vuonghm/.cshrc.common",__LINE__);
my $out_fn;my %input_elements; my $webaddon="-u www-data ";
die "$usage\n" if (!defined $opt_f);
if (defined $opt_o){
	$out_fn=$opt_o;
}else{
	$out_fn="$opt_f.config.ini";
}
my $date=`date "+%Y%m%d"`;chomp $date;
$webaddon='' if (!$opt_w);
#################Variable Declarations and zip formats ####################
my $BASE_DIR;
my $web_bin=`dirname $0`;chomp $web_bin;
my $annovar_dir='/bioinfoC/hue/annovar/current/';
my $mirna_exe='/bioinfoC/hue/mirna/mirna_impact_all.pl';
################Read in user inputs########################
my $web_input=$opt_f;my $web_dir=`dirname $opt_f`;chomp $web_dir;
if (!$web_dir){
	if ($bin=~/(avia1|dev)/){
		$web_dir='/mnt/webrepo/fr-s-abcc-avia1/public/complete';
	}else{
		$web_dir='/mnt/webrepo/fr-s-abcc-avia0/public/complete';
	}
}
##This is where we append user config files for the web server
my $user_config=`grep ^user.email= $web_input`;chomp $user_config;
my $email=($user_config=~/email=(\S+)\@/)?$1:'';
sendMsg ("ERROR","no email could be found $user_config\n",'') if (!$email);
if ( `grep modulesid $web_input`=~/[0-46]$/ && $email ){
	my $id=$1;my $addon='';
	if ($id==6){
		$addon="_cascade";
	}
	if ($email && -e "$web_dir/userspecific/$email/config.ini$addon"){
		SAFE_PROCESS ("cp $web_input $web_input.orig; cat $web_dir/userspecific/$email/config.ini$addon >> $web_input\n",__LINE__);
	}
}
my ($isFatal_msg,$err_msg,$isLocal);$isLocal=0;
SAFE_PROCESS("$perl $web_bin/convertTxtToConfig.pl -f $web_input -o $out_fn",__LINE__);#generates $out_fn
################read in config ######################
my %cfg;
tie %cfg , 'Config::IniFiles', (-file =>"$out_fn");
my $scripts_bin=$cfg{'Utilities'}{'scripts_dir'};chomp $scripts_bin;
my $workdir="$cfg{'Target_Info'}{'base_dir'}";chomp $workdir;
my $original_dir=`dirname $0`;chomp $original_dir;
$original_dir.="/../";
if (exists $cfg{'UserInfo'}{'ipaddress'} && $cfg{'UserInfo'}{'ipaddress'} ne ''){
}
SAFE_PROCESS("mkdir $original_dir/data/$cfg{'Target_Info'}{'label'};chmod 777 $original_dir/data/$cfg{'Target_Info'}{'label'}",__LINE__) if (! -e "$original_dir/data/$cfg{'Target_Info'}{'label'}");
my $af=0;#calculate allele frequency from vcf file for jason lih
################## Move to AVA work directory #############################
SAFE_PROCESS ("mkdir $workdir \n",__LINE__) if (! -e "$workdir");
open (LOG,">$workdir/web_wrapper.log") or sendMsg("ERROR, write err", "cannot open $workdir/webwrapper.log for writing\n",'');
open (STDERR, ">>&LOG"); open (STDOUT,">>&LOG");
print STDERR "Redirecting stderr and stdout to $workdir/web_wrapper.log\n";
if ((-e "$workdir/$out_fn") && (`diff $out_fn $workdir/$out_fn | wc -l` >1)){#session id will always be different
	my $ct=1;my $curr="$workdir/$out_fn";
	until (! -e "$workdir/$out_fn\.$ct") {
		if (`diff $out_fn $workdir/$out_fn.$ct | wc -l `>0){
			$ct++;
		}else{
			last;
		}
	}
	SAFE_PROCESS("cp $curr $workdir/$out_fn\.$ct\n",$0,__LINE__);
}
SAFE_PROCESS ("mv $out_fn $workdir/config.ini\n",__LINE__);
################## Launch AVA #############################
chdir ("$workdir");
my $program=$cfg{'UserInfo'}{'program'};
my $runcmd=printCmdFile();
print "About to run $runcmd \n";
my $gid=SAFE_PROCESS("cd $workdir;$runcmd",__LINE__);
if (!$gid){
	sendMsg( "ERROR, QSUB submission", "could not submit the $runcmd to the grid for processing\nPlease check director $workdir at $0 LN".__LINE__."\n","");
	exit;
}
my $isDone="";		my $convert=`which convert`;chomp $convert;
my $date=`date "+%Y%m%d"`;chomp $date;

#if [ -s traces ]
if ( -e "$workdir/traces" && `ls $workdir/traces | wc -l`>0){
#then
#isDone=''
	$isDone="";
#$convert $workdir/*png $cfg{'Target_Info'}{'label'}.pdf
	SAFE_PROCESS("convert $workdir/traces/*png $workdir/traces/$cfg{'Target_Info'}{'label'}.pdf", __LINE__);
#cp $workdir/traces/$cfg{'Target_Info'}{'label'}.pdf cfg{'Target_Info'}{'base_dir'}/traces/$cfg{'Target_Info'}{'label'}.csv $original_dir/data/$cfg{'Target_Info'}{'label'}
	SAFE_PROCESS("mv $workdir/AVA_$cfg{'Target_Info'}{'gene'}\_$date.csv $workdir/$cfg{'Target_Info'}{'label'}.csv",__LINE__);
	SAFE_PROCESS("cp $workdir/traces/$cfg{'Target_Info'}{'label'}.pdf $workdir/$cfg{'Target_Info'}{'label'}.csv $original_dir/data/$cfg{'Target_Info'}{'label'}", __LINE__);
}elsif ($program=~/ftp_upload|trace_archives_download/){
	$isDone=" -z ";
	SAFE_PROCESS("cp $workdir/$cfg{'Target_Info'}{'label'}.csv $original_dir/data/$cfg{'Target_Info'}{'label'}", __LINE__);
}elsif ($program=~/impact_analysis_only/i){
#	my $final_output_name;
#	if (-e "$workdir/annovar_wrpr.output" && ! -z "$workdir/annovar_wrpr.output"){
#		if ($isLocal){
#			if (exists $cfg{'UserInfo'}{'out_name'}){
#				$final_output_name=$cfg{'UserInfo'}{'out_name'};chomp $final_output_name;
#			}else{
#				$final_output_name="annovar_wrpr.output";
#			}
#			my $filelocation=`dirname  $cfg{'UserInfo'}{'input_fullpath'}`;chomp $filelocation;
#			SAFE_PROCESS("mv $workdir/annovar_wrpr.output $filelocation/$final_output_name\n",__LINE__);
#			SAFE_PROCESS("ln -s $filelocation/$final_output_name $workdir/annovar_wrpr.output\n",__LINE__);
#			SAFE_PROCESS("ln -s $filelocation/$final_output_name $original_dir/data/$cfg{'Target_Info'}{'label'}/$cfg{'Target_Info'}{'label'}.tsv \n",__LINE__);
#			SAFE_PROCESS("cp $workdir/config.ini $original_dir/data/$cfg{'Target_Info'}{'label'}/. \n",__LINE__);
#		}else{
#			SAFE_PROCESS("mv $workdir/annovar_wrpr.output $original_dir/data/$cfg{'Target_Info'}{'label'}/$cfg{'Target_Info'}{'label'}.tsv\n",__LINE__);
#			SAFE_PROCESS("ln -s $original_dir/data/$cfg{'Target_Info'}{'label'}/$cfg{'Target_Info'}{'label'}.tsv annovar_wrpr.output\n",__LINE__);#make a link in the original directory
#			SAFE_PROCESS("cp $workdir/config.ini $original_dir/data/$cfg{'Target_Info'}{'label'}/. \n",__LINE__);
#		}
#		SAFE_PROCESS("touch $workdir/DONE.processing\n",__LINE__);
#	}else{
#		$isDone=" -z ";
#	}
#	sendMsg("Thank you for using the AVIA software" , "Your analysis $cfg{'Target_Info'}{'label'} is now complete.  You can directly link to your page by clicking below or by cutting the link below and pasting into any web browser:\n".
#		"http://avia.abcc.ncifcrf.gov/apps/site/results/?id=$cfg{'Target_Info'}{'label'}\n\nYou can also retrieve any of your other submissions by using our data retrieval page at : ".
#		"http://avia.abcc.ncifcrf.gov/apps/site/retrieve_a_request to see your results and input your id and analysis type.  Your results will be stored for 1 week from the date of submission","$email");
#	SAFE_PROCESS ("echo '$cfg{'Target_Info'}{'label'}'  >> $original_dir/update.db\n");
}elsif ($program=~/(viz|ia_upload)/){
	print "exiting";
	exit;
}else{
	die "what's this program at ln".__LINE__."\n";
}
#SAFE_PROCESS("$perl $web_bin/uploadPublicToDB.pl -c $workdir/config.ini $isDone" , __LINE__);
#SAFE_PROCESS("touch $workdir/DONE.dbupdate",__LINE__);
exit(0);
sub printCmdFile{#no longer accepts LMT inputs
#BTW you are already in the uniqid folder on tork!!!! /bioinfoC/AVA/PUBLIC_SUBMISSIONS/uniqid!!!!!
	my $cmd;$refver='';
	my $_user_original;
	if ($program=~/(ftp_upload|trace_archives_download)/i){
		open (BAT, ">$workdir/runAVIA.bat" ) || sendMsg ("ERROR, File I/O", "Cannot open runAVIA.bat in $workdir\n","");
		print BAT "#!/bin/bash\nPBS_O_WORKDIR=$workdir\n". 
			"cd \${PBS_O_WORKDIR}\numask 0000\n";
		#use search specifications to download chromats :)
		my $cmd=qq($perl $scripts_bin/get_traces_db.pl -c $workdir/config.ini -g
printenv  >> retrieval.log
touch DONE.processing
	);
		if ($program=~/trace_archives_download/i){
			if (exists $cfg{'UserInfo'}{'trace_reseq_only'}){
				$cmd="$perl $scripts_bin/prepareTraceArcQuery.pl -c config.ini -r -l 100000\n".$cmd;
			}else{
				$cmd="$perl $scripts_bin/prepareTraceArcQuery.pl -c config.ini -l 10000\n".$cmd;
			}
		}else{
			SAFE_PROCESS ("cp $cfg{'UserInfo'}{'file'} $workdir/.",__LINE__);
			my $filename=`basename $cfg{'UserInfo'}{'file'}`;chomp $filename;
			unPack("$workdir/$filename");
		}
			print BAT $cmd;
			close BAT;
			SAFE_PROCESS ("chmod +x $workdir/runAVIA.bat\n",__LINE__);
			return ("$qsub $webaddon -l pvmem=500mb -e $workdir/stderr -o $workdir/stdout -q grande -M vuonghm\@mail.nih.gov $workdir/runAVIA.bat");
	}elsif ($program=~/(impact_analysis_only|viz|mirna)/i){
		my ($filelocation,$filename);
		open (BAT, ">$workdir/runAVIA.bat" ) || sendMsg ("ERROR, File I/O", "Cannot open runAVIA.bat in $workdir\n","");
		print BAT "#!/bin/bash\nPBS_O_WORKDIR=$workdir\n". 
			"cd \${PBS_O_WORKDIR}\numask 0000\n";
		if (exists $cfg{'UserInfo'}{'file'}){
			#this is for uploaded files; need to move off of avia webserver as SGE cannot access these directories
			$filelocation=`dirname $cfg{'UserInfo'}{'file'}`;
			$filename=`basename $cfg{'UserInfo'}{'file'}`;
			$isLocal=0;
		}elsif (exists $cfg{'UserInfo'}{'input_fullpath'}){
			#this is for internal users; already tork accessible so they just need a link
			$filelocation=`dirname  $cfg{'UserInfo'}{'input_fullpath'}`;
			$filename=`basename  $cfg{'UserInfo'}{'input_fullpath'}`;
			$isLocal=1;
		}elsif (exists $cfg{'UserInfo'}{'input_path'}){
			die "Haven't coded yet\n";
		}else{
			sendMsg("ERROR", "Could not find file or input_fullpath in your config file\n");
		}
		chomp $filelocation;
		chomp $filename;
		$_user_original=$filename;$_user_original=~s/^(\d{14}\.)//;
		#check that the original file is not zipped or tarred
		if ($filename=~/(zip|tar|gz)$/){
			#untar change filename and change the config file, but do not rename the file
			$_rename=0;
			$filename=pop(@{checkAndSanitize($cfg{'UserInfo'}{'file'})});
		}
		if ($isLocal){
			my $fileflag=0;
			open (TEST,"<$filelocation/$filename") or $fileflag=1;
			close TEST if (!$fileflag);
			if (!-e "$filelocation/$filename" || $fileflag){
				sendMsg("ERROR, Non-existent input file", 
				"You have submitted an incorrect path and file to our AVIA server.  Your request cannot be processed ($cfg{'Target_Info'}{'label'})\n",$email);
			} 
			SAFE_PROCESS("ln -s $filelocation/$filename .\n",__LINE__);
			
		}else{
			#this is on the server and 
			if (exists $cfg{'UserInfo'}{'wastyped'} && $cfg{'UserInfo'}{'wastyped'}){
				
				if (`grep -P '\t' $filelocation/$filename | wc -l`==0){
					open (FILE,"<$filelocation/$filename" ) or die "Cannot open $filelocation/$filename\n";
					open (OUT,">$workdir/$filename" ) or die "Cannot open $workdir/$filename\n";
					while (my $line=<FILE>){
						$line=~s/\s{1,}/\t/g;
						print OUT "$line\n";
					}
					close OUT;close FILE;
				}else{
					SAFE_PROCESS ("/usr/bin/dos2unix $filelocation/$filename\n") if (-e "/usr/bin/dos2unix");
					SAFE_PROCESS( "cp $filelocation/$filename $workdir/$filename \n",__LINE__);
				}
			}else{
				SAFE_PROCESS( "cp $filelocation/$filename $workdir/$filename \n",__LINE__);
			}
		}
		if ($program=~/mirna/){
			$cmd="$perl $mirna_exe $workdir $filename\n";
			$web_input=`basename $web_input`;chomp $web_input;
			$cmd.="if [ -e \"$workdir/Impact_Table\" ]\nthen\ncp Impact_Table $original_dir/data/$cfg{'Target_Info'}{'label'}/$cfg{'Target_Info'}{'label'}\n$perl $bin/sendmsg.pl $original_dir/complete/$web_input 1\nelse\n$perl $bin/sendmsg.pl $original_dir/complete/$web_input 0\nfi\n";
			print BAT "$cmd\n";
			close BAT;
			SAFE_PROCESS ("chmod +x $workdir/runAVIA.bat\n",__LINE__);
			return ("$qsub $webaddon -l pvmem=1GB -e $workdir/stderr -o $workdir/stdout -q small -M vuonghm\@mail.nih.gov $workdir/runAVIA.bat \n");
		}else{#this is so we can run the annovar script on the grid and so that the next scripts do not begin to run until it is done
			my $run=qq(starttime=`date +%s`
run=`echo \$runid|sed 's/.abcc1.ncifcrf.gov//g'`
diff=0
myrun=`qstat | grep \$run|wc -l`;
until  [ \$myrun = 0  ] ;do
        sleep 30;
        currtime=`date +%s`;
        diff=\$(( (\$currtime-\$starttime ) ));
        if [ \$diff -gt 3600 ];then
                break
        fi
        myrun=`qstat | grep \$run|wc -l`;
done);
			my $type=$cfg{'UserInfo'}{'inputformat'};my $annotated=0;
			my $convert_exe="$perl $annovar_dir/convert2annovar.pl -format ";
			if ($type=~/(pileup|cg|gff3-solid)/i){
				#need to convert to annovar format;
				SAFE_PROCESS ("$convert_exe $1  -outfile ANNOVAR.input $filename \n",__LINE__);
				$filename="ANNOVAR.input";
			}elsif ($type=~/(soap|maq|vcf4|casava|clcbio|modhgvs|hgvs)/i){
				SAFE_PROCESS ("$convert_exe $1  --includeinfo --outfile ANNOVAR.input $filename\n",__LINE__);
				$filename="ANNOVAR.input";
			}elsif ($type=~/(vcf|annovar)/){
				#annot_vcf,vcf,annot_anvr,anvr
				if ($type=~/^vcf$/){
					my $addon= "-allallele ";
					if (exists $cfg{'UserInfo'}{'calcAF'} && $cfg{'UserInfo'}{'calcAF'}=~/true/i){
						$addon.=" -af ";
						$af=1;
					}else{
						$addon='';
					}
					SAFE_PROCESS ("$convert_exe vcf4  --includeinfo $addon --outfile ANNOVAR.input $filename\n",__LINE__);
					#check the output file for af
					if ($af && -e "cvrt2anvr.stderr.log"){
						my $pass=` grep -e no_af -e afct cvrt2anvr.stderr.log | wc -l`;
						$af=1 if ($pass>0);
					}
					$filename="ANNOVAR.input";
				}
				if ($type=~/annot/){
					$annotated=1;
				}
			}
			#make a file with the database to run against
			my $db_loc="/SeqIdx/annovardb/";
			if ($cfg{'UserInfo'}{'ver'}=~/mm/){
				$db_loc.="mousedb/";$dbinfofile=~s/XXXX/mousedb/;
				$refver=$cfg{'UserInfo'}{'ver'};
			}else{
				$db_loc.="humandb/";$dbinfofile=~s/XXXX/humandb/;
				
				if ($cfg{'UserInfo'}{'ver'}=~/hg18/){
					#convert here
					print STDERR "[INFO] Running /bioinfoB/dwnld/ucsc/tools/liftOver $filename /bioinfoC/bobs/alt_splice_shalini/ctcf/hg18ToHg19.over.chain hg19_converted.tmp hg19_failed -positions...\n";
					SAFE_PROCESS ("/bioinfoB/dwnld/ucsc/tools/liftOver $filename /bioinfoC/bobs/alt_splice_shalini/ctcf/hg18ToHg19.over.chain hg19_converted.tmp hg19_failed -positions\n");
					my $newfilename=`find . -maxdepth 1 |grep bedmapped\$ `;chomp $newfilename;
					if (-e $newfilename && -z $newfilename){
						sendMsg("ERROR converting hg18 to hg19","$filename was improperly formatted or has no valid conversions\n",'');
						exit;
					}
					SAFE_PROCESS( "mv $filename $filename.orig\n");
					SAFE_PROCESS ("mv $newfilename hg19_ANNOVAR.input\n");
					system ("echo '#$cfg{'UserInfo'}{'ver'}_origPos' > hg19_hg18converted.txt\n");
					print STDERR "[INFO] Running: cut -f1-3 $filename.org | sed 's/[\\t]/:/' | sed 's/[\\t]/-/' >> hg19_hg18converted.txt\n";
					SAFE_PROCESS (" cut -f1-3 $filename.orig | sed 's/[\\t]/:/' | sed 's/[\\t]/-/' >> hg19_hg18converted.txt\n");
					SAFE_PROCESS ("rm liftOver* hg19_converted.tmp -rf\n");
					$filename='hg19_ANNOVAR.input';
					$refver='hg19';
				}else{
					$refver='hg19';
				}
			}
			open (DBFILE,">$workdir/searchTheseDBs.txt" ) ;#$workdir/
			$_rename=1;#we want to rename the databases in the format <ORG>_NAME.TXT for ANNOVAR
			if ($af){
				print DBFILE "$workdir/ANNOVAR.input.hg19_af\n";
			}
			if ($cfg{'UserInfo'}{'ver'} ne "$refver"){
				print DBFILE "$workdir/hg19_hg18converted.txt\n";
			}
			DBFILE: foreach my $dbfile (keys %{$cfg{'UserInfo'}}){
				next if ($dbfile !~/annotdb_(.*)/ && $dbfile!~/^userdefined_db/);#db_snp132=on
				my $db=$1;
				if ($db ne '' ){#these are the standard databases available in annovar; they are in the /SeqIdx/annovardb/<org>db/ directories
					if ($db!~/userdefined/i){
						push(@files, "$db_loc/$refver\_$db.txt") ;#let annovar_wrpr take care of the dbs that do not exist for the specified organism
					}else{
						next;
					}
				}else{#these are the user uploaded files #userdefined_db1=myfile.txt
					my $file=`basename $cfg{'UserInfo'}{$dbfile}`;chomp $file;
					#move to the work directory that is accessible on the grid and so we can uncompress if necessary
					SAFE_PROCESS ("/usr/bin/dos2unix $cfg{'UserInfo'}{$dbfile}\n");
					SAFE_PROCESS  ("cp $cfg{'UserInfo'}{$dbfile} $workdir/$file\n",__LINE__);#this is the zip/tar files, better to move around compressed files
					#if the file is not a zip/tar format
					if ($file!~/(zip|tar|tgz|gz)/){
						if (! -e "$file" || `file $file`!~/ASCII.*text/i){
							#do not process anything that isn't a text file
							next;
						}
						$file=abcc_rename($file);
						push(@files, "$workdir/$file");
					}else{
						#move the file to the final directory
						my $arr=checkAndSanitize($file);#an array of filenames is returned in case there are more than one files
						if ($#{$arr}==-1){#no acceptable files could be extracted
							sendMsg ("ERROR:checkAndSanitize", "Could not execute checkAndSanitize on $file\n$?");die;
						}else{
							foreach my $em (@{$arr}){
								push(@files, "$workdir/$em");#see next note for why we need to include the \n here
							}
							###do not change the next line...logic foreach of the tar'd files, we add it to the list of databases to check.  
							###But these should not have any scores because there is no way to get the user score filter to this 
							###line at this time, score is requred to have the format: \d+,\d+ and not just one column for annovar
							next DBFILE;
						}
					}
				}
				if (($cfg{'UserInfo'}{'filter'}=~/yes/i) && $dbfile=~/(sift|polyphen)/i ){
					$files[$#files].="\t".$cfg{'UserInfo'}{"$1\_cutoff"};
				}
			}
			reorder(\*DBFILE);
			my $filter_addon='';
			if ($cfg{'UserInfo'}{'modulesid'}==6){
				#do some filtering first;
				open (FILORDER,">$workdir/filter_order.txt") or die "Cannot open filter_order.txt";
				foreach my $ids (sort keys %{$cfg{'UserInfo'}}){
					if ($ids=~/^filter(\d+)/){
						print FILORDER "$cfg{'UserInfo'}{$ids}\t". $cfg{'UserInfo'}{"keep$1"};
						if (exists $cfg{'UserInfo'}{"keep$1cutoff"}){
							print FILORDER "\t".$cfg{'UserInfo'}{"keep$1cutoff"}."\n";
						}else{
							print FILORDER "\t0\n";
						}
						
					}
					
				}
				close FILORDER;
				$filter_addon=" -F filter_order.txt ";
			}
			my $ensembl='';
			if ($cfg{"UserInfo"}{'useEnsembl'}=~/on/){
				$ensembl= "-E";
			}
			$cmd= "runid=`$perl $annovar_dir/annovar_qsub_wrpr.pl $ensembl -i $filename $filter_addon -f searchTheseDBs.txt -b -o $refver -W`\nBAZ\n";my $othercmd;
			my $runfiltergenelist=getFromCfg(\%cfg);
			if ($program=~/viz/){
				my $nextFile;
				#do some other stuff and add to cmd
				my $success=addViz(\*DBFILE, $filename);
				my $addon;
				if ($cfg{'UserInfo'}{'strict'}=~/yes/){
					$addon.=" -x ";
				}
				
				if  ( $success && -f "$workdir/$success" ){
					$filename=$success;
					$cmd= "runid=`$perl $annovar_dir/annovar_qsub_wrpr.pl $ensembl -i $filename $filter_addon -f searchTheseDBs.txt -b -o $refver $addon -W`\nBAZ\n";
					$nextFile="$filename.annovar_wrpr.output";
				}elsif ($annotated){
					$cmd="";#use that file as the next step of the analysis
					$nextFile=$filename;
				}else{
					die "You were unsuccessful in parsing data for subtractive analysis at line ($success)".__LINE__."\n";
				}
				my $sfx='';
				if ($runfiltergenelist){
					$cmd.="\n$perl $annovar_dir/filtergenelist.pl -f genelist -c $out_fn -i $nextFile";
					if ($cfg{'UserInfo'}{'genelist_type'}=~/filter/i){
						$sfx=".orig";#we do this because of the highlight vs filter option
					}
					$cmd.="\ncp $nextFile.orig genelist $original_dir/data/$cfg{'Target_Info'}{'label'}/.";#so this will be in the zip file
				}
				my $runcircos=0;
				foreach my $checker (keys %{$cfg{'UserInfo'}}){
					if ($cfg{'UserInfo'}{$checker}=~/(^on$|^true$)/i){
						$runcircos++;#print "working on $checker,$cfg{'UserInfo'}{$checker}...pass1";<STDIN>;
					}elsif ($checker=~/^key_/ && $cfg{'UserInfo'}{$checker}=~/(pop|all)/){#precomputed
						$runcircos++;#print "working on $checker,$cfg{'UserInfo'}{$checker}...pass2";<STDIN>;
					}
					last if $runcircos;
				}
				if ($runcircos ){
					$cmd.="\n$perl /users/abcc/vuonghm/scripts.dir/parseData/run_viz_pipeline.pl -f $nextFile$sfx -c config.ini -d $workdir -s -N ";
					$othercmd="cp config.ini viz/*png viz/circos.input viz/forthelegend.txt $original_dir/data/$cfg{'Target_Info'}{'label'}/."; 
					$othercmd.="\nif [ -e \"viz/subtractive.out\" ]\nthen\ncp viz/subtractive.out $original_dir/data/$cfg{'Target_Info'}{'label'}/$cfg{'Target_Info'}{'label'}_genicOnly.tsv\nfi";
				}elsif ( $cfg{'UserInfo'}{'modulesid'}==2 || $cfg{'UserInfo'}{'modulesid'}==1){
					$cmd.="\n$perl /bioinfoC/hue/scripts.dir/parseData/run_viz_pipeline.pl -F $nextFile$sfx -c config.ini -d $workdir -s -N ";
					$othercmd="cp config.ini viz/*png viz/circos.input viz/forthelegend.txt $original_dir/data/$cfg{'Target_Info'}{'label'}/."; 
					$othercmd.="\nif [ -e \"viz/subtractive.out\" ]\nthen\ncp viz/subtractive.out $original_dir/data/$cfg{'Target_Info'}{'label'}/$cfg{'Target_Info'}{'label'}_genicOnly.tsv\nfi";
				}
				$cmd.="\n$perl $annovar_dir/rename_headers.pl $nextFile";
			}elsif ($program=~/impact_analysis_only/){
				my $nextFile="$filename.annovar_wrpr.output";
				$cmd.="\n$perl /users/abcc/vuonghm/scripts.dir/parseData/run_viz_pipeline.pl -f $nextFile -c config.ini -d $workdir -s -N ";
			}
			close DBFILE;
			
			$cmd=~s/\-b\s\-o/\-b \-z \-o/ if (exists $cfg{'UserInfo'}{'show_flanking'} );
			$cmd.="\n$perl /bioinfoC/hue/scripts.dir/parseData/summarize.pl $filename $_user_original";
			$cmd.="\ncp $filename.annovar_wrpr.output $original_dir/data/$cfg{'Target_Info'}{'label'}/$cfg{'Target_Info'}{'label'}.tsv;$othercmd";
			my $user_orig_fn=$filename;
			$user_orig_fn=~s/(^\d{14}\.|\.txt$|\.csv$|\.tsv$|\.vcf$)//g;#if this was a duplicate file and I added the timestamp to the beginning
			if ($cfg{'UserInfo'}{'name_column_header'}==1){
				$cmd.="\nif [ -e \".$filename.annovar_wrpr.output\" ]\nthen\ncp .$filename.annovar_wrpr.output $original_dir/data/$cfg{'Target_Info'}{'label'}/$cfg{'Target_Info'}{'label'}.tsv\nfi";
			}
			$cmd.="\nif [ -e \"$filename.STATS\" ]\nthen\ncp $filename.STATS $original_dir/data/$cfg{'Target_Info'}{'label'}/$cfg{'Target_Info'}{'label'}.stats\nfi";
			foreach my $otherfiles ("sift.txt","pp2.txt", "DD.txt","$filename.jpg","$filename.STATS.jpg","hg19_failed", "avia.database.info.txt"){
				my $newname;
				if ($otherfiles=~/$filename(.*)$/){
					$newname=$user_orig_fn.$1;
				}else{
					$newname=".";
				}
				if ($otherfiles=~/avia.database/){#this is empty but we need this in order for the web script running.php to recognize that the process is completed
					$cmd.="\nif [ ! -e \"$otherfiles\" ]\nthen\ntouch $otherfiles\nfi";
				}
				$cmd.="\nif [ -e \"$otherfiles\" ]\nthen\ncp $otherfiles $original_dir/data/$cfg{'Target_Info'}{'label'}/$newname\nfi";
			}
			$cmd=~s/BAZ/$run/;
			print BAT "$cmd\nchmod 777 * -R\nchmod 777 $original_dir/data/$cfg{'Target_Info'}{'label'} -R";
			close BAT;
			SAFE_PROCESS ("chmod +x $workdir/runAVIA.bat\n",__LINE__);
			return ("$qsub -l pvmem=1GB -e $workdir/stderr $webaddon -o $workdir/stdout -q small -m ba -M vuonghm\@mail.nih.gov $workdir/runAVIA.bat \n");
		}#ends if/else $program=~/mirna/
	}else{
		die "Cannot distinguish the program type ($program)\n";
	}
	return  ("$qsub $webaddon -l pvmem=500mb -e $workdir/stderr -o $workdir/stdout -q small -m ba -M vuonghm\@mail.nih.gov $workdir/runAVIA.bat\n");
}
sub addViz{
	#first check subtractive inputs
	my $fh=shift;my $fn=shift;
	my $href=$cfg{'UserInfo'};#shortcut;
	my $targetfile=$fn;
	my $addl_annot=0;
	if (exists $$href{'cgi_pop'} && $$href{'cgi_pop'}){
		my @arr=split(",",$$href{'cgi_pop'});
		foreach my $pop (@arr){
			$addl_annot++;
			print $fh "/SeqIdx/annovardb/humandb/hg19_$pop.txt\n";
		}
	}elsif (exists $$href{'groupfile'} && $$href{'groupfile'} ){
		#user has populations data in their input file
		system ("cp $$href{'groupfile'} .\n");
		my $opts;
		if ($$href{'inputformat'}=~/vcf/){
			$opts="-f $fn ";
		}else{
			$opts="-F $fn "
		}
		#for now we assume 2 populations in the data file
		SAFE_PROCESS("$perl /bioinfoC/hue/scripts.dir/parseData/parsePopulationData.pl $opts -t specific -p ". `basename $$href{'groupfile'}`, __LINE__);
		#should produce 2 files hg19_<pop1>_hom_only.txt and hg19_<pop2>_hom_only.txt
		print "LN".__LINE__."cut -f2 $$href{'groupfile'} |grep -ve $$href{'target_pop'} |head -n1\n";
		my $other=`cut -f2 $$href{'groupfile'} |grep -ve $$href{'target_pop'} |head -n1`;chomp $other;
		$targetfile = "hg19_$$href{'target_pop'}\_hom_only.txt";
		my $otherfile=`ls *$other\_hom_only.txt `;chomp $otherfile;
		print $fh "$workdir/$otherfile\n";
		$addl_annot++;
	}
	if ($addl_annot){
		return $targetfile;
	}else{
		return $fn;
	}
}
sub SAFE_PROCESS {
	my $cmd=shift;
	my $line_nbr=shift;
	chomp $cmd if ($cmd=~/\n$/);
	if ($cmd=~/qsub/){
		my $id=`$cmd`;chomp $id;
		if ($id=~/(\d+.abcc1)/){
			return $1;
		}elsif ($id=~/(\d+)/){
			print LOG "Grid submission entry: $id\n";
			return $1;
		}else{
			sendMsg("ERROR, QSUB", "execution of command from $0 ($cfg{'Target_Info'}{'label'})at line $line_nbr FAILED\n$?");
			die;
		}
	}else{
		eval {
			#print STDERR "$cmd\n";
			system ("$cmd 2>&1\n");	
		};
		if ($?){
			sendMsg ("ERROR, PROCESSING", "execution of $cmd from $0 ($cfg{'Target_Info'}{'label'}) at line $line_nbr FAILED\n$?");
			die;
		}else{
			return 1;
		}
	}
}
sub sendMsg{
	my $subject=shift;
	my $msg=shift;
	my $email_addy=shift;
	if ($email_addy ne '' ){
		$ADMIN.=",$email_addy";
		
	}
	my  %mail = ( To  =>   "$ADMIN",
			BCC => 'vuonghm@mail.nih.gov',
            From    => 'AVIA_admin@mail.nih.gov',
            Subject=> "$subject",
            Message => "\n$msg\n"
           );

  sendmail(%mail) or die $Mail::Sendmail::error;
  if ($subject=~/error/i){
  		print STDERR "$msg\n";
  		exit (1);
  }else{
		return;
  }
}
sub abcc_rename{
	my $prename=shift;
	if (!$_rename){return $prename;}
	my $orig=$prename;
	if ($prename!~/^($cfg{'UserInfo'}{'ver'}|$refver)\_/){
		$prename=$cfg{'UserInfo'}{'ver'}."_$prename";
	}
	if ($prename!~/.txt/){
		$prename.=".txt";
	}
	if ($prename ne "$orig"){
		system ("mv $orig $prename\n");
	}
	return $prename;
}
sub checkAndSanitize{#NOTE: we use this for both databases and the user input file
	my $xfile=shift;
	my $xcmd='';my @xarr;
	if ($xfile=~/bz2{0,1}$/){
		unPack($xfile);
		$xfile=~s/.bz2{0,1}$//g;
		push(@xarr,abcc_rename($xfile));
	}elsif ($xfile=~/(.(tar\.|t){0}gz)$/){
		print "Found $1|->\t$2|->\t$3|\n";
		my $sufx=$1;
		unPack($xfile);
		$xfile=~s/$sufx//;
		push(@xarr,abcc_rename($xfile));
	}elsif ($xfile=~/zip$/){
		#find all files in the zipfile
		my @features=split("\n",`$unzip -l $xfile`);
		unPack($xfile);
		for (my $i=3;$i<=$#features-2;$i++){
			print "$i...working on features $features[$i]....\n";
			if ($features[$i]=~/(\S+)$/){
				my $old=$1;
				my $newname=abcc_rename($old);
				system ("mv $workdir/$old $workdir/$newname\n") if ($old ne "$newname");
				push (@xarr,$newname);
			}else{
				die "wtf $features[$i]...\n";
			}
		}
	}elsif ($xfile=~/tar$/){
		print "$xfile\n\t$tar -tvf $xfile\n";
		foreach my $file (split("\n",`$tar -tvf $xfile`)){
			chomp $file;
			if ($file=~/(\S+)$/){
				next if ($file=~/\/$/);#directory
				$file=$1;
				my $info=`$tar -xvf $xfile --get $file`;
				$info=`file $info`;chomp $info;
				if ($info=~/ASCII.*text/i){
					print STDERR "[Info] ok:mv'ing $file .\n";
					system ("mv $file .\n");#so what if it's the same file.  This is to get rid of any additional directories in the archive
					push (@xarr,abcc_rename($file));
				}else{
					system ("rm $file -fr\n");
				}
			}
		}
	}elsif($xfile=~/(tar.|t)gz$/){
		print STDERR "haven't coded for $1 yet at $0 LN".__LINE__."\n";
	}else{
		print STDERR "Haven't coded for this suffix...skipping\n";
	}
	return (\@xarr);
	
}
sub unPack{
	my $xfile=shift;my $newdir;
	sendMsg("ERROR, UserInput FileName", "This ($xfile) is going to cause problems on unix, rename and resubmit\n", "") if ($xfile=~/\s/);
	if($xfile=~/^\//){
		$newdir=`dirname $xfile`;chomp $newdir;
		chdir $newdir;
	}
	if ($xfile=~/(ab1|scf|trace|txt)$/){
		return 1;
	}elsif ($xfile=~/tar.gz$/ || $xfile=~/tgz$/ ){
		system ("$gunzip < $xfile | $tar xvf -\n");
	}elsif ($xfile=~/tar$/){
		system ("$tar -xf $xfile\n");
	}elsif ($xfile=~/bzip$/){
		system ("$bzip $xfile\n");		
	}elsif ($xfile=~/(gz|z)$/i){
		system ("$gunzip $xfile\n");
	}elsif ($xfile=~/zip$/i){
		system ("$unzip $xfile\n");
	}elsif ($xfile eq ''){
		die "haven't coded yet\n";
		#need to unpack by going into each of the directories and unpacking one at a time
	}else{
		sendMsg( "ERROR, USER INPUT extension", "The user has specified a unknown package ($xfile)...please code", "");
	}
	chdir("$workdir");#changed this from original directory on 12/12/11
	return;
}

sub monitor {
	my $grid_id=shift;
	my $start_time=`date "+%s"`;chomp $start_time;my $curr_time=$start_time;
	my $MONITOR_LEN=600;$start_time+=$MONITOR_LEN;#roughly 1 hour
	my $notDone=1;
	unless ( $start_time<=$curr_time ){
		$notDone=`qstat | grep $grid_id`;chomp $notDone;
		if ($notDone){
			$curr_time=`date "+%s"`;
			sleep (int($MONITOR_LEN/10));
		}
	}
	return $notDone;
}

sub reorder{
	my $fh=shift;
	my $flag=0;
	print LOG "In Reorder...\n";
	if (! -e $dbinfofile){
		print LOG "$dbinfofile does not exist\n";
		my $mysql_exe=`which mysql`;chomp $mysql_exe;
		if ( $mysql_exe != ''){
			$dbinfofile="$web_dir/database.info";
			system ("$mysql_exe -h sqldb1.abcc -u abccruser -pabccrpwd <$web_dir/sql >$dbinfofile\n");
			open (DBDUMP,"<$dbinfofile") or return;
		}else{
			return;
		}
	}
	open (DBDUMP,"<$dbinfofile") or $flag=1;my $max=-1;my %priOrder;
	print LOG "Successfully opened a database file to rearrange ($dbinfofile)\n";
	my $nidx=my $priority=-1;
	while (<DBDUMP>){
		my @elements=split("\t",$_);chomp $elements[$#elements];
		if ($_=~/^database\_id/){
			for (my $i=0;$i<=$#elements;$i++){
				if ($elements[$i]=~/database\_category\_idx/){
					$priority=$i;
					last;
				}elsif ($elements[$i]=~/database\_annvr\_name/){
					$nidx=$i;
				}
			}
		}else{
			die "headers not found" if ($nidx<0 || $priority<0);
			$priOrder{$elements[$nidx]}=$elements[$priority];
			$max=$elements[$priority] if ($max<$elements[$priority]);
		}
	}
	close DBDUMP;
	$max++;
	my %ordered;my $chkcount=$#files;
	print Dumper (\%priOrder);
	foreach my $dbname (@files){
		if ($dbname=~/$refver\_(.*)\.txt/ ){
			my $db=$1;
			if (exists $priOrder{$db}){
				$ordered{$priOrder{$1}}.="$dbname\n";
			}else{
				print "$db not exist in priOrder hash!\n";
				$ordered{$max}.="$dbname\n";
			}
		}else{
			print "$opt_o\n";
			$ordered{$max}.="$dbname\n";
		}
	}
	print Dumper (\%ordered);
	@files=();
	foreach my $key (sort {$a<=>$b} keys %ordered){
		print $fh "$ordered{$key}";
	}
}