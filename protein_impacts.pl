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
my $MAX_PROCESSES=4;
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
		-p path to your databases
		-o Output name DEFAULT \$opt_I/i.hg19_PTM.txt
		-h <if specified, no header>
		-D Full Path to Database
);
$opt_v||="hg19";
my $org=($opt_v=~/hg\d+/)?"humandb":"mousedb";
my $noheader=($opt_h)?1:0;
die $usage if (!defined $opt_i && !defined $opt_I);
my $input;
if ($opt_i){
	$input=$opt_i;
}elsif($opt_I){
	$input=$opt_I;
}
if (!$opt_o){
	my $tmp=$input;
	$tmp=~s/(exonic\_){0,}variant\_function(.collapsed){0,}//;
	$opt_o="$tmp$opt_v\_PTM.txt";
}
my $keepline=($opt_k)?1:0;
my $database_directories=(defined $opt_p)?$opt_p:"/SeqIdx/annovardb/$org/";
#the list of files to be queried; path relative to $database_directories
my @files2Query;
my $database;
my %prot_db;
my $output_header;
if (!$opt_D){
	 %prot_db=(
		'PhosphoSitePlus.regulatorySites'=>{ 
				'file'=>'hg19_PhosphoSitePlus.regulatorySites.txt',
				'href'=>1,
			},
		'PhosphoSitePlus.SumoylationSites'=>{ 
				'file'=>'hg19_PhosphoSitePlus_SumoylationSites.txt',
				'href'=>1,
			},
		'PhosphoSitePlus.MethylationSites'=>{ 
				'file'=>'hg19_PhosphoSitePlus_MethylationSites.txt',
				'href'=>1,
			},
		'PhosphoSitePlus.O-GlcNAcSites'=>{ 
				'file'=>'hg19_PhosphoSitePlus_O-GlcNAcSites.txt',
				'href'=>1,
			},
		'PhosphoSitePlus.AcetylationSites'=>{ 
				'file'=>'hg19_PhosphoSitePlus_AcetylationSites.txt',
				'href'=>1,
			},
		'PhosphoSitePlus.UbiquitinationSites'=>{ 
				'file'=>'hg19_PhosphoSitePlus_UbiquitinationSites.txt',
				'href'=>1,
			},
		'PhosphoSitePlus.DiseaseSites'=>{ 
				'file'=>'hg19_PhosphoSitePlus_Disease-associatedSites.txt',
				'href'=>1,
			},
		'PhosphoSitePlus.KinaseSubstrateSites'=>{ 
				'file'=>'hg19_PhosphoSitePlus_Kinase-SubstrateSites.txt',
				'href'=>1,
			},
		'PhosphoSitePlus.PhosphorylationSites'=>{ 
				'file'=>'hg19_PhosphoSitePlus_PhosphorylationSites.txt',
				'href'=>1,
			},
		'Phosida'=>{
				'file'=>'hg19_phosida_sites.txt',
				'href'=>1,
			}
	);
	foreach my $prot_fn (keys %prot_db){
		open (REG,"$database_directories/prot/$prot_db{$prot_fn}{'file'}") or die "Cannot open $prot_db{$prot_fn}{'file'} file\n";
		my %href=();
		while (<REG>){
			my ($pp,$info)=split("\t",$_);chomp $info;
			$href{$pp}=$info;
		}	
		$prot_db{$prot_fn}{'href'}=\%href;
		close REG;
	}
	$output_header="#PTM";
}else{
	if ($opt_D!~/^\//){
		$opt_D=$database_directories."prot/$opt_D";
	}
	$output_header=`basename $opt_D`;chomp $output_header;
	$output_header=~s/($opt_v\_|\.txt.*$)//g;
	open (REG,"<$opt_D") or die "Cannot open $opt_D file\n";
	my %href=();
	while (<REG>){
		my ($pp,$info)=split("\t",$_);chomp $info;
		$href{$pp}=$info;
	}	
	$prot_db{$opt_D}{'href'}=\%href;
	close REG;
}	
#print Dumper (\%prot_db);die;
die "$database_directories does not exist...please check and try again or specify opt_d\n" if (! -e "$database_directories");
my %report;
if ($opt_i){
	open (INPUT,"<$opt_i") or die "Cannot open $opt_i\n";
}elsif ($opt_I){
	open (INPUT,"<$opt_I") or die "Cannot open $opt_i\n";
}
print "done opening input\n";
#my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
#$pm->start and next; # do the fork
#$pm->wait_all_children;
open (OUT,">$opt_o") or die "Cannot open $opt_o for writing\n";
my $index= my $feat_idx=-1;my $output="";my $linenbr=0;
LINE: while (my $line=<INPUT>){
	chomp $line;
	if ($line=~/unknown/){
		print OUT "-\n" if ($keepline);next;
	}
	my ($type,$info,@other,@elements);
	if ($opt_i){
		($linenbr,$type,$info,@other)=split("\t",$line);
	}else{
		$linenbr++;
		@elements=split("\t",$line);
		if ($line=~/^(.*Rel\spos\sto\sfeature)/){
			my $foo=$1;
			while ($foo=~/\t/g){$index++;}
			$index++;
			$feat_idx=$index-2;
			if ($keepline && !$noheader){print OUT "#$output_header\n";}
			next LINE;
		}elsif ($index<0){
			if ($keepline && $linenbr==1 && !$noheader){print OUT "#PTM\n";}
			$index=2;$feat_idx=0;
		}else{
			print OUT "$output" and $output='' if ($linenbr%100000==0);
			if ($elements[$feat_idx]!~/exonic/ || $elements[$feat_idx]=~/(\*|ncRNA)/i || $elements[$index]=~/unknown/i){#asterisks denotes that info was not found in the exonic_variants_function
				$output.="-\n" if ($keepline) ;
				next LINE;
			} 	
		}
		$info=$elements[$index];
	}
	my (@trx_info)=split(",",$info);
	my ($gene,$target);
	my $found='';
	foreach my $trx (@trx_info){
		next if (!$trx);
		my @elements=split(":",$trx);
		if (!$gene){
			if ($elements[0]=~/(synon|frameshift|stop)/){
				$gene=$elements[1];
			}else{
				$gene=$elements[0];
			}
		}
		if ($trx=~/p\.(\w\d+)([A-Z]{1,2})/i){
			$target=uc("$gene:$1");my $protChange=$2;
			if ($target=~/:$protChange/){#synonymous
				next;
			}
			foreach my $database (sort keys %prot_db){
				my $reg=$prot_db{$database}{'href'};
				$found.=(exists $$reg{$target} && $found!~/$$reg{$target}/)?"$$reg{$target},":"";
			}
		}elsif ($trx=~/p\.(\d+)\_(\d+)\w+:(\w{3})/){
			my ($start,$stop,$triplet,$moon)=($1,$2,$3,$4);
			my $aa=$cdn2prot{uc($triplet)};
			$target=uc("$gene:$aa$start");
			foreach my $database (sort keys %prot_db){
				my $reg=$prot_db{$database}{'href'};
				$found.=(exists $$reg{$target} && $found!~/$$reg{$target}/)?"$$reg{$target},":"";
			}
		}elsif ($trx=~/(wholegene|NR\_)/){
		}else{
			warn "couldn't parse protein info from line $linenbr($trx)$elements[$feat_idx]\t$line";
		}
	}
	if ($found){
		chop $found if ($found=~/,$/);
		$output.= "$found\n" ;
	}elsif ($keepline){
		$output.="-\n";
	}
	
}
print OUT $output if ($output);
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
