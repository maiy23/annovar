#!/usr/bin/perl
=head
This script 
=cut
use strict;
use FindBin;
use lib "$FindBin::Bin/";
use Data::Dumper;
use Getopt::Std;
use vars qw($opt_c $opt_h $opt_e $opt_l $opt_w);
getopts("c:e:l:wh");
umask(0000);
my %cfg;
my $usage=qq(
	$0 
	REQUIRED:
	-c <config file> 
	OPTIONAL:
	[-e <email address of user>] 
	[-l <AVIA id>]
	[-w If specified than it is a web user and we need to transfer to web directory]

	Your configuration file looks like: Optional values are in brackets [ ]

	INPUT_TYPE=<vcf4|bed|anvr|clcbio|varscan2>
	INPUT_DIRECTORY=<full path to the directory with all your input files>
	[DATABASE_FILE=<full path to the database file> ] #This is your searchDBs.txt file;default is set
	[SUFFIX]=<txt|vcf> ##default ending
	[LOG]=<fullpath to name of log file>

);

print $usage and exit if ($opt_h);
###These are standard executables
my $annovar_dir=`dirname $0`;chomp $annovar_dir;#changed for FDI'/bioinfoC/hue/annovar/prod/';##"prod" directory will not change ; "current" directory will have the most recent AVIA run
my $annovar_exe="$annovar_dir".'/annovar_qsub_wrpr.pl';
my $annovar_convert="$annovar_dir/".'convert2annovar.pl';
##Check input variables
if (!defined $opt_c){
	print "You must specify a configuration file\n$usage\n";
	exit;
}
##read in configuration file
open (CONFIG,"<$opt_c") or die "Cannot open config file($opt_c)\n";
while (my $line=<CONFIG>){
	$line=~s/\s//g;
	my ($key,$value)=split("=",$line);chomp $value;
	$cfg{uc($key)}=$value;
}
$cfg{'DATABASE_FILE'}||='/bioinfoC/AVA/FDI/searchDBs.txt';
my $dbfn=`basename $cfg{DATABASE_FILE}`;chomp $dbfn;
close CONFIG;

die "One or more of your inputs are not valid\n$usage" if (!exists $cfg{'INPUT_TYPE'} || ! exists $cfg{'INPUT_DIRECTORY'} || !-e $cfg{'INPUT_DIRECTORY'});
my $logfile=(exists $cfg{LOG})?$cfg{LOG}:"$cfg{INPUT_DIRECTORY}/log.err";
if ($logfile!~/^\//){
	$logfile=$cfg{INPUT_DIRECTORY}."/$logfile";
}
open (LOG,">$logfile") or die "Cannot open log file $logfile\n";
open (STDERR, '>&LOG') or die "Cannot redirect STDERR to log file $logfile\n";
open (STDOUT, '>&LOG') or die "Cannot redirect STDOUT to log file $logfile\n";
print Dumper (\%cfg);
##start evaluating directories
chdir($cfg{'INPUT_DIRECTORY'});
my $sfx=".txt";
# read directory
opendir(my $DIR, $cfg{INPUT_DIRECTORY}) or die "Error opening $cfg{INPUT_DIRECTORY}\n";
# Loop for each file in the directory
my %files;
while (my $file = readdir($DIR)){
	next if ($file=~/^(\.*|PROC|$opt_c|$dbfn|log.err|.*bat|.*log|config.ini)$/);
	next if (! -f $file);
	if (!exists $files{$file}){
		$files{$file}=$file;
	   $files{$file}=~s/\.(vcf|txt)$//;
	   $sfx=$1;
	   if (! -e "$cfg{INPUT_DIRECTORY}/PROC/$files{$file}"){
	   		SAFE_PROCESS  ("mkdir -p PROC/$files{$file}\n");
	   		chdir ("PROC/$files{$file}");
	   		system ("ln -s ../../$file $file\n");#link the new file
	   		system ("ln -s $cfg{DATABASE_FILE} .\n");#link the list of databases file
	   		SAFE_PROCESS ("perl $annovar_convert -format $cfg{INPUT_TYPE} -o ANNOVAR.input -allallele --includeinfo $file\n",__LINE__);#run converter
	   		chdir ($cfg{INPUT_DIRECTORY});
	   }else{
	   		print STDERR "[INFO] PROC/$files{$file} already exists! Skipping conversion and copy step...\n";
	   }
	}
}
closedir($DIR);
print Dumper (\%files);
# run annovar's convert script
chdir("PROC");
system ("mkdir All \n") if (!-e 'All');
SAFE_PROCESS ("cut -f1-5 */ANNOVAR.input |sort -u >> All/All.input.txt\n",__LINE__);
chdir ("All");
system ("ln -s $cfg{DATABASE_FILE} .\n");
my $whattomonitor=monitor ("perl $annovar_exe -f $dbfn -i All.input.txt",0);#the second argument 0 tells the subroutine to return and not wait for the job to complete; returns the grid job id
##make a global bat file in case the grid job failed....
open (BAT,">post_processing.bat") or die "Cannot open post_processing.bat\n";
print BAT qq(
#!/bin/sh
PBS_O_WORKDIR=$cfg{INPUT_DIRECTORY}/PROC
cd \$PBS_O_WORKDIR
#PBS -l cput=290:000:00,pcput=290:00:00,walltime=290:00:00
#PBS -j oe -l pvmem=2GB,mem=2GB,ncpus=2 -e \$PBS_O_WORKDIR/PROC/stderr -m e
perl $annovar_dir/mergeOneAnnot_multiInputs.pl $cfg{INPUT_DIRECTORY}/PROC/All/All.input.txt.annovar_wrpr.output $cfg{INPUT_DIRECTORY}/PROC
);
if ($cfg{INPUT_TYPE}=~/(vcf|varscan)/i){
print BAT qq(for i in `ls |grep -ve '^All\$'`
do
	cd \$i
	ln -s ../All/avia.database.info.txt .
	perl $annovar_dir/writeVCF.pl -f vcf -w -o \${i}_withannot.vcf -I ANNOVAR.input.annovar_wrpr.output -p \${i}.$sfx -A $dbfn
	cp \${i}_withannot.vcf $cfg{INPUT_DIRECTORY}
	sfx='withannot.vcf'
	chmod 777 $cfg{INPUT_DIRECTORY}/. -R
	cd ../
done
);
}else{
	print BAT qq(for i in `ls |grep -ve '^All\$'`
do
	cd \$i
	cp ANNOVAR.input.annovar_wrpr.output $cfg{INPUT_DIRECTORY}/\$i.withannot.output
	sfx='withannot.output'
	cd ..
done
if [ -e "PROC/All/avia.database.info.txt"] ;then
	cp PROC/All/avia.database.info.txt .
	chmod 777 $cfg{INPUT_DIRECTORY}/* -R
fi
);
}
if (defined $opt_e && !$opt_w){
	my $label=(defined $opt_l)?$opt_l:$$;
print BAT qq(if [ -e '$cfg{INPUT_DIRECTORY}/PROC/All/DONE.All.input.txt.anvr.wrpr' ];then
	perl $annovar_dir/sendmsg.pl '$opt_e' '$label completed processing' 'Thank you for using AVIA.  AVIA has completed processing your internal submission.  Please view your annotation files in $cfg{INPUT_DIRECTORY}.  If you do not have access to this directory, please contact the AVIA team at NCI-FrederickAVIA\@mail.nih.gov.'
else 
	perl $annovar_dir/sendmsg.pl '$opt_e' '$label errored out' 'Thank you for using AVIA.  AVIA has errored out during processing your internal submission.  Please view our FAQs for internal submissions to debug what went wrong at /bioinfoA/scripts/avia/README.  Or please contact the AVIA team at NCI-FrederickAVIA\@mail.nih.gov.'
fi
);
}elsif (defined $opt_e && $opt_w){
	my $label=$opt_l;
	#Send users a different message and move stuff to the web directory
	my ($dev,$server);
	if ($opt_l=~/dev/){
		$server='/mnt/webrepo/fr-s-abcc-avia0/public/data';$dev='dev';
	}else{
		$server='/mnt/webrepo/fr-s-abcc-avia0/public/data';$dev='';
	}
	print BAT qq(if [ -e '$cfg{INPUT_DIRECTORY}/PROC/All/DONE.All.input.txt.anvr.wrpr' ];then
	cd $cfg{INPUT_DIRECTORY};zip $label.zip *\${sfx}
	cp searchTheseDBs.txt $label.zip  avia.database.info.txt $server/$label
	perl $annovar_dir/sendmsg.pl '$opt_e' '$label completed processing' 'Thank you for using AVIA.  AVIA has completed processing your submission.  Please view your annotation files at:

http://avia$dev.abcc.ncifcrf.gov/apps/site/download/?file=$label/$label.zip
'
else 
	perl $annovar_dir/sendmsg.pl 'vuonghm\@mail.nih.gov' '$label errored out' 'AVIA has errored out during processing your internal submission.  Please view our FAQs for internal submissions to debug what went wrong at /bioinfoA/scripts/avia/README.  Or please contact the AVIA team at NCI-FrederickAVIA\@mail.nih.gov.'
fi		
		);
}
close BAT;
system ("chmod 777 post_processing.bat\n");
# aggregate data across all samples ## this will greatly speed up the amount of grid time
# segregate data by sample after ANNOVAR has run
#convert back to VCF if VCF4 was original input
my $success=0;
if ($whattomonitor){
	$success=monitor($whattomonitor,1);
}
if ($success){
	monitor("qsub $cfg{INPUT_DIRECTORY}/PROC/All/post_processing.bat\n",0);
}else{
	print "Your job id $whattomonitor did not complete or there was an error.  When resolved, please run post_processing.bat\n";
}

close LOG;
sub SAFE_PROCESS{
	my $cmd=shift;
	my $line=shift;
	eval{
		print STDERR "[INFO] Running $cmd\n";
		system ($cmd);
	};
	if ($?){
		die "Could not run $cmd at line $line\n";
	}
	return;
}
sub monitor{#we do this because the grid directive to run another file does not work!!
	my $targets=shift;
	my $monitor=shift;
	my $done=1;
	if ($monitor){
		my $timeout=`date "+%d%H%M%S"`;chomp $timeout;$timeout+=3600000;#1hour
		print STDERR "Monitoring $targets....\n";
		until (`$targets` eq '' || !$done){
			sleep (60);
			my $curr_time=`date "+%d%H%M%S"`;chomp $curr_time;
			if ( $curr_time>$timeout){$done=0;}
		};
		return $done;
	}else{
		my $gid=`$targets`;
		if ($gid=~/^(\d{6,})/){
			return "qstat | grep $1";
		}else{
			return "$gid\n";
		}
	}
}