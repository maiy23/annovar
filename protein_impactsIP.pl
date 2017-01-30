	#!/usr/bin/perl
=head
Reads the annovar exonic_variantsfunction file and queries protein coordinates files
=cut
use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/Parallel-ForkManager-0.7.9/lib";
use Parallel::ForkManager;
my (%aas,%cdn2prot);
my $MAX_PROCESSES=5;
my %reg;
use vars qw( $opt_p $opt_o $opt_i $opt_d $opt_I $opt_v $opt_k $opt_h $opt_D);
getopts("i:o:d:p:I:v:f:D:kh");
my $usage=qq($0 
	[ONE OF THE FOLLOWING MUST BE SPECIFIED]
		-i <INPUT FILE> *exonic_variant_function from ANNOVAR
		-I <input file> annovar_wrapper file  (collapsed)
	[OPTIONAL]
		-k keep all lines in order of input (for -I option only);
		-v organism DEFAULT hg19; 
		-p path to your databases #default /SeqIdx/annovardb/humandb/prot
			Note:  all databases listed must be in the same directory, if not, this script must be run separately and this option changed
		-o Output name DEFAULT \$opt_I/i.hg19_XXX
		-h <if specified, no header>
		-D Databases to run #default only PTM is run
);
$opt_v||="hg19";
my $org=($opt_v=~/hg\d+/)?"humandb":"mousedb";
my $noheader=($opt_h)?1:0;
die $usage if (!defined $opt_I);
my $input;
if ($opt_i){
	die "This is deprecated!\n";exit;
}elsif($opt_I){
	$input=$opt_I;
}
# We use the XXX to change the name of each of the databases
if (!$opt_o){
	my $tmp=$input;
	$tmp=~s/(exonic\_){0,}variant\_function(.collapsed){0,}//;
	$opt_o="$tmp$opt_v\_XXX";
}else{
	$opt_o=~s/.txt/XXX/;
	if ($opt_o!~/XXX/){
		$opt_o.="XXX";
	}
}
my $keepline=($opt_k)?1:0;
my $database_directories=(defined $opt_p)?$opt_p:"/SeqIdx/annovardb/$org/prot/";
#the list of files to be queried; path relative to $database_directories
my @files2Query;
my %prot_db;
my $output_header;
my @databases;
if ($opt_D){
	@databases=split(",",$opt_D);
}else{
	push(@databases,"PTM");
}
die "$database_directories does not exist...please check and try again or specify opt_d\n"  if (! -e "$database_directories");
foreach my $database (@databases){
	if ( $opt_v=~/hg19/ && $database=~/(\bPTM\b|PhosphoSitePlus)/ ){
		$database=$1;
		 my %ptm_prot_db=(
			'PhosphoSitePlus.regulatorySites'=>{ 
					'file'=>'hg19_PhosphoSitePlus.regulatorySites.txt',
				},
			'PhosphoSitePlus.SumoylationSites'=>{ 
					'file'=>'hg19_PhosphoSitePlus_SumoylationSites.txt',
				},
			'PhosphoSitePlus.MethylationSites'=>{ 
					'file'=>'hg19_PhosphoSitePlus_MethylationSites.txt',
				},
			'PhosphoSitePlus.O-GlcNAcSites'=>{ 
					'file'=>'hg19_PhosphoSitePlus_O-GlcNAcSites.txt',
				},
			'PhosphoSitePlus.AcetylationSites'=>{ 
					'file'=>'hg19_PhosphoSitePlus_AcetylationSites.txt',
				},
			'PhosphoSitePlus.UbiquitinationSites'=>{ 
					'file'=>'hg19_PhosphoSitePlus_UbiquitinationSites.txt',
				},
			'PhosphoSitePlus.DiseaseSites'=>{ 
					'file'=>'hg19_PhosphoSitePlus_Disease-associatedSites.txt',
				},
			'PhosphoSitePlus.KinaseSubstrateSites'=>{ 
					'file'=>'hg19_PhosphoSitePlus_Kinase-SubstrateSites.txt',
				},
			'PhosphoSitePlus.PhosphorylationSites'=>{ 
					'file'=>'hg19_PhosphoSitePlus_PhosphorylationSites.txt',
				},
			'Phosida'=>{
					'file'=>'hg19_phosida_sites.txt',
				}
		);
		my %href=();
		foreach my $prot_fn (keys %ptm_prot_db){
			# print "working on $prot_fn($database_directories/$prot_db{$prot_fn}{'file'})\n";
			open (REG,"$database_directories/$ptm_prot_db{$prot_fn}{'file'}") or die "Cannot open $ptm_prot_db{$prot_fn}{'file'} file at LINE ". __LINE__."\n";
			while (<REG>){
				my ($pp,$info)=split("\t",$_);chomp $info;
				if (exists $href{$pp}){
					$href{$pp}.=";" if ($href{$pp}!~/;$/);
					$href{$pp}.="$info";
				}else{
					$href{$pp}=$info;
				}
			}	
			close REG;
		}
		$prot_db{$database}=\%href;
	}else{
		my $name="$database_directories/$opt_v\_$database.txt";
		print "[DEBUG] opening $name\n";
		if (!-e "$name"){
			print "$name does not exist...try $database_directories/$database...\n";
			$name="$database_directories/$database";#user specified and may not conform to AVIA
		}
		open (REG,"<$name") or die "Cannot open $name file at LINE ". __LINE__."\n";
		my %href=();
		while (<REG>){
			my ($pp,$info)=split("\t",$_);chomp $info;
			$href{uc($pp)}=$info;
		}	
		$prot_db{$database}=\%href;
		close REG;
	}	
}
my %report;
open (INPUT,"<$opt_I") or die "Cannot open $opt_i\n";
print "done! Now opening input\n";
#my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
#$pm->start and next; # do the fork
#$pm->wait_all_children;
my $index= my $feat_idx=-1;my $output="";my $linenbr=-1;
my %found;my %blankline;
LINE: while (my $line=<INPUT>){
	chomp $line;#print "working on $line\n";<STDIN>;
	my ($type,$info,@other,@elements);
	
	$linenbr++;
	@elements=split("\t",$line);
	if ($line=~/^(.*Rel\spos\sto\sfeature)/){
		my $foo=$1;
		while ($foo=~/\t/g){$index++;}
		$index++;
		$feat_idx=$index-2;
		# if ($keepline && !$noheader){$blankline{$linenbr}=1;}
		next LINE;
	}elsif ($index<0){
		die "Cannot have index <1";
		# if ($keepline && $linenbr==1 && !$noheader){$found{$linenbr}='-';}
		$index=2;$feat_idx=0;
	}elsif ($elements[$feat_idx]!~/exonic/ || $elements[$feat_idx]=~/(\*|ncRNA)/i || $elements[$index]=~/unknown/i){#asterisks denotes that info was not found in the exonic_variants_function
		next LINE;
	}
	$info=$elements[$index];
	my (@trx_info)=split(",",$info);
	my ($gene,$target);
	foreach my $trx (@trx_info){
		next if (!$trx);
		my @elements=split(":",$trx);
		if (!$gene){
			if ($elements[0]=~/(synon|frameshift|stop|start)/){
				$gene=$elements[1];
			}else{
				$gene=$elements[0];
			}
		}
		if ($trx=~/p\.(\w\d+)([A-Z]{1,2})/i){
			$target=uc("$gene:$1");my $protChange=$2;
			foreach my $database (keys %prot_db){
				 if (exists $prot_db{$database}{$target}){
				 	$found{$linenbr}{$database}.="$prot_db{$database}{$target}," if ($found{$linenbr}{$database}!~/$prot_db{$database}{"$target$protChange"},/)
				 }elsif(exists $prot_db{$database}{"$target$protChange"} ){
				 		$found{$linenbr}{$database}.=$prot_db{$database}{"$target$protChange"}."," if ($found{$linenbr}{$database}!~/$prot_db{$database}{"$target$protChange"},/);
				 }
			}
		}elsif ($trx=~/p\.(\d+)\_(\d+)\w+:(\w{3})/){
			my ($start,$stop,$triplet,$moon)=($1,$2,$3,$4);
			my $aa=$cdn2prot{uc($triplet)};
			$target=uc("$gene:$aa$start");#print "working with " .__LINE__.":$target($trx)\n";
			foreach my $database (keys %prot_db){
				$found{$linenbr}{$database}.="$prot_db{$database}{$target}," if (exists $prot_db{$database}{$target});
			}
		}elsif ($trx=~/(wholegene|NR\_)/){
		}else{
			# warn "couldn't parse protein info from line $linenbr($trx)$elements[$feat_idx]\t$line";
		}
	}	
}
# $output_header=`basename $database`;chomp $output_header;
# 		$output_header=~s/($opt_v\_|\.txt.*$)//g;
# print Dumper (\%found);
my @fh;
# print $linenbr;
foreach my $fh (keys %prot_db){
	my $fn=$opt_o;$fn=~s/XXX/$fh/;
	open (FH,">$fn") or die "Cannot open $fn\n";
	if (!$noheader){print FH "#$fh\n";}#0 is header
	for (my $i=1;$i<=$linenbr;$i++){
		if (!exists $found{$i}){
			print FH "-\n";
		}elsif (exists $found{$i}{$fh} && $found{$i}{$fh}=~/\w/){
			chop $found{$i}{$fh} if ($found{$i}{$fh}=~/,$/);
			print FH "$found{$i}{$fh}\n";
		}else{
			print FH "-\n";
		}
	}
}
close INPUT;
close OUT;
sub _defineAAs {
    %aas=(
    "A","Ala",    	"R","Arg",
    "N","Asn",		"D","Asp",
    "C","Cys",		"Q","Gln",
    "E","Glu",		"G","Gly",
    "H","His",		"I","Ile",
    "START","start","L","Leu",
    "K","Lys",		"M","Met",
    "F","Phe",		"P","Pro",
    "S","Ser",		"T","Thr",
    "W","Trp",		"Y","Tyr",
    "V","Val",		"*","stop",
    "","-", "STOP","stop","X","stop",
    "-","frameshift");
	 %cdn2prot=("AAA","K","AAC","N","AAG","K","AAT","N",                                            
          "ACA","T","ACC","T","ACG","T","ACT","T",                                            
          "AGA","R","AGC","S","AGG","R","AGT","S",                                            
          "ATA","I","ATC","I","ATG","M","ATT","I",                                            
          "CAA","Q","CAC","H","CAG","Q","CAT","H",                                            
          "CCA","P","CCC","P","CCG","P","CCT","P",                                            
          "CGA","R","CGC","R","CGG","R","CGT","R",                                            
          "CTA","L","CTC","L","CTG","L","CTT","L",                                            
          "GAA","E","GAC","D","GAG","E","GAT","D",                                            
          "GCA","A","GCC","A","GCG","A","GCT","A",                                            
          "GGA","G","GGC","G","GGG","G","GGT","G",                                            
          "GTA","V","GTC","V","GTG","V","GTT","V",                                            
          "TAA","*","TAC","Y","TAG","*","TAT","Y",                                            
          "TCA","S","TCC","S","TCG","S","TCT","S",                                            
          "TGA","*","TGC","C","TGG","W","TGT","C",                                            
          "TTA","L","TTC","F","TTG","L","TTT","F");                                           
    return;
}

sub cvrtCDN2AA{
	my $triplet=shift;
	if (exists $aas{$cdn2prot{$triplet}}){
		return $aas{$cdn2prot{$triplet}};
	}else{
		die "could not convert triplet $triplet\n";
	}
}
