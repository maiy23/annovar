#!/usr/local/bin/perl

=head
v1.0 Original 
This helper script was designed to convert original format back to vcf
v1.5 
-Use of collapsed ANNOVAR inputs i.e. one file instead of variant_function and exonic_variant function
v2.0
-Different algorithm
-added the -a option, which specifies that the -allalleles option was used during conversion of vcf->annovar format
=cut
use strict;
use Switch;
use Getopt::Std;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use lib "$FindBin::Bin/Parallel-ForkManager-0.7.9/lib";
use Parallel::ForkManager;
use ABCC;
use vars qw($opt_O $opt_F $opt_i $opt_o $opt_t $opt_f $opt_p $opt_a $opt_n $opt_b $opt_A $opt_I $opt_w);
getopts("i:t:f:t:p:a:n:A:o:I:O:F:bw");
my $usage=qq(
	$0 -p <original file>
	One of the following must be specified:
		-i <ANNOVAR input>  #prefix for ANNOVAR.input.variant_function and ANNOVAR.input.exonic_variant_function
		-I <collapsed ANNOVAR input>  ##annovar_wrpr.txt
	OPTIONAL	
		-a <additional database names> || -A <filename with list of database names>
			-comma separated
		-c <columns for annot>  ##not yet implemented
			Used if you are not using ANNOVAR
		-f <original format of file>
			This was the original file format of the input file before convert2annovar.pl
		-n <name column-use this file as the ID in the ID column>	
			only used with opt_t:rsid and append option
			Note: specify opt_b if the name file DOES NOT HAVE A HEADER!!! (not off by one);
		-o <output filename>
			DEFAULT append with .vcf
		-t <user_defined|rsid|append|none|>
			This is to determine what to write into the id column of the vcf file \(col 3\)
			Can be used in conjunction with -n option \(above\)
		-w <with abcc_tag tags> extra processing on data integrity
			specify this option this if you used --allallele option when running convert2annovar.pl
	
	Note: <ANNOVAR input> is the basename for the input for ANNOVAR.  We will use the
	ANNOVAR.variant_function and ANNOVAR.exonic_variant_function files as annotation.  If
	you have an annotation file from another source, then it must match the original file
	in line number and you will have to specify which columns are used to annotate.
);
my $debug=0;
my ($annovar_base,$annot_fn,$OffbyOne,$FunSeq2,%fdi_headers,$exac,%exac);
if ( (!defined $opt_p && $opt_t=~/vcf/i) || (!defined ($opt_i) && !defined $opt_I && !defined $opt_F)){ die $usage;}
$opt_t ||= 'user_defined';$opt_t=lc($opt_t);
$opt_f ||= 'vcf4';
$opt_O ||= 'hg19';
if ($opt_b){
	$OffbyOne=0;
}else{
	$OffbyOne=1;
}
if (defined $opt_i){
	if ($opt_i=~/.collapsed/){
		#do not run collapse script;
		$annovar_base=$opt_i;
		$annot_fn=$opt_i;
		$annovar_base=~s/.variant_function.collapsed//g;
	}elsif ($opt_i=~/.variant_function/){
		die "One or more of your ANNOVAR output files do not exist!\n" if (!-e $opt_i || !-e "$opt_i.exonic_variant_function");
		collapseANNOVAR($opt_i);
		$annovar_base=~s/.variant_function//g;
		$annot_fn=$opt_i.".collapsed";
	}
}elsif (defined $opt_I){
	$annot_fn=$opt_I;
	$annovar_base=$opt_I;
}elsif (defined $opt_F){
	$annot_fn=$opt_F;
	$annovar_base="$opt_F.$opt_O"."_";
}else{
	die "One or more of your ANNOVAR output files do not exist!\n" if (!-e "$opt_i.variant_function" || !-e "$opt_i.exonic_variant_function");
	collapseANNOVAR("$opt_i.variant_function");
	$annovar_base=$opt_i;
	$annot_fn=$opt_i.".variant_function.collapsed";
}
$opt_o ||= "$annovar_base.vcf";
my $FILE_LOAD=10000;my  %hash;my ($annot_href,$name_href,$headers_ref,@FunSeq2_headers);my $chunk;
################################################
my @db_searched_files;
if (defined $opt_a){
	@db_searched_files=split(",",$opt_a);
}elsif (defined $opt_A){
	open (DBS,"<$opt_A" ) || die "Cannot open $opt_A\n";
	while (my $file=<DBS>){
		chomp $file;
		$file=`basename $file`;chomp $file;$file=~s/.txt//g;
		$file="$annovar_base.$file";
		if (-e $file && !-z $file){
			push(@db_searched_files,$file);
		}elsif ($opt_I || $opt_F){
			push(@db_searched_files,$file);
		}else{
			#ignoring becase it failed ANNOVAR.
		}
	}
	close DBS;
}
if ($opt_t=~/(rsid|append)/ && !defined $opt_n){
	my $err;
	if (defined $opt_a && $opt_a=~/snp/i && $opt_a=~/,{0,1}([^,]*),{0,1}/i  ){
		$opt_n=$1;
		if ($opt_n=~/^\s*$/){$err=1;}else{
		print STDERR "You are using $opt_n as your identifier\n";}
	}else{$err=1;
	}
	if ($err){
		print STDERR "[WARN] Could not find your SNP file to use as ID in the vcf output...using info in the original as ID\n";
		undef($opt_n);$opt_t='user_defined';
	}
}
if (defined $opt_n && $opt_a!~/$opt_n/i){
	push (@db_searched_files,$opt_n);print "pushing $opt_n\n";
}
print "opening output file ($opt_o) for writing\n";
open (FILE,">$opt_o") or die "Cannot open $opt_o\n";
	#check that the original file and the transformed file have the same number of elements
	#this is no longer a valid check to tri-allelic inputs in vcf file
#print "[INFO] Checking input files for compatibility (this may take a while)\n";
#my $orig_ct=`grep -ve '#' $opt_p | wc -l`;chomp $orig_ct;
#my $transform_ct=`wc -l $annot_fn`;chomp $transform_ct;$transform_ct=~s/($annot_fn|\s)//g;
#die "New file:$transform_ct and the original:$orig_ct do not have the same number of lines...Please double check\n" if ($transform_ct!=($orig_ct+1));
#print STDERR "[Info] Both input files match lines($transform_ct).\n";
my $multiAllele=0;#This keeps track of whether the =(allele,allele2) needs to be printed for a particular line (only if there is a multi allele present on line and $opt_w is set)
open (ORIG,"<$opt_p") or die "cannot open ($opt_p) for reading\n";
open (ANNOT, "<$annot_fn") or die "$! $annot_fn\n";
print "opening $annot_fn for reading...\n";
my %idx;my $ct=0;my $found=0;
  while (!eof(ANNOT) and !eof(ORIG)) {
    my $orig = <ORIG>;
  	 if ($orig=~/^#CHROM/ && $opt_f=~/vcf/i){
  	 	if (!$found ){
			my $headers=findDBheaders();
			print FILE $headers;
			$found=1;
		}
		print FILE $orig;
	}elsif ($orig=~/^#/ ){
		if ($opt_f=~/vcf/i){
			if ($orig=~/##INFO=/ && !$found ){##incomplete header and not true VCF but we have to put this somewhere!
				#here to enter more header info from your searched database;
				my $headers=findDBheaders();
				print FILE $headers;
				$found=1;
				# print Dumper (\%fdi_headers);<STDIN>;
			}
			print FILE $orig;
		}else{
			my ($header,$footer)=getAuxHeaders();
			#printing VCF header...
			print FILE "$header". findDBheaders()."$footer\n";#this does not have a new line at the end
			$ct++;
			
		}
	}else{
		my $annot=<ANNOT>;
		$annot=~s/[\n\r]//g;
		if ( ($opt_F && $annot=~/(Summary\t){0,}Variant ID/) || ($annot=~/^#/)){
#			if (-z $opt_o){
#				my ($header,$footer)=getAuxHeaders();
#				#printing VCF header...
#				print FILE "$header". findDBheaders()."$footer\n";#this does not have a new line at the end
#				$ct++;
#			}
			my ($href,$aref)=readAnnotBlockFDI($annot_fn,1,1);
			for (my $i=0;$i<=$#{$aref};$i++){
				next if ($$aref[$i] eq '');
				$$aref[$i]=~s/#//g;
				if ($$aref[$i]=~/^ANNOVAR annot/i){
					$idx{'ANVR_TYPE'}=$i;
				}elsif ($FunSeq2 && $$aref[$i]=~/ENCODE.ANNOTATED;/i){
					$$aref[$i]=~s/#//g;
					@FunSeq2_headers=split(";",$$aref[$i]);
					$idx{$$aref[$i]}=$i;
				}elsif ($$aref[$i]=~/^Gene$/i){
					$idx{'ANVR_GENE'}=$i;
				}elsif ($$aref[$i]=~/^Gene[\s|_]ori$/i){
					$idx{'GENE_ORI'}=$i;
				}elsif ($$aref[$i]=~/(Rel pos|Annot Feat)/i){
					$idx{'ANVR_ANNOT'}=$i;
				}elsif ($$aref[$i]=~/Comment/i){
					$idx{'last_col'}=$i;
		#					$i++ if ($i<$#{$aref});#the next column is the column with sample information and not just a comment
				}elsif (exists $fdi_headers{$$aref[$i]}){
					# print 'LN'.__LINE__." mapped: $fdi_headers{$$aref[$i]}\n";
					$idx{$fdi_headers{$$aref[$i]}}=$i;
				}elsif ($$aref[$i]=~/#{0,1}(.*)/i ){
					my $tmp=$1;
					next if $tmp=~/^(chrom|comment|samplename|prot.pos)$/i;
					# print 'LN'.__LINE__."$1\n";#print Dumper (\%fdi_headers);
					$idx{uc($tmp)}=$i;
				}else{
					print "Skipping $$aref[$i]\n";
				}
			}
			$annot=<ANNOT>;
			$annot=~s/[\n\r]//g;
			# print Dumper (\%fdi_headers);
			# print Dumper (\%idx);<STDIN>;
		}
		my @orig_info=split("\t",$orig);
		if ($opt_f=~/vcf/i){
			my @multi_alleles;$multiAllele=0;
			if ($orig_info[4]=~/,/i && $opt_f=~/vcf/i && $opt_w){
				@multi_alleles=split(",",$orig_info[4]);
				$multiAllele=1;
			}elsif (!$opt_w && $opt_f=~/vcf/i &&  $orig_info[4]=~/^([\w\-]+),/i  ){
				push(@multi_alleles,$1);
			}else{
				push (@multi_alleles,$orig_info[4]);
			}
			my $aref=\@orig_info;
			my $last_annot;
			ALLELE: for (my $i=0;$i<=$#multi_alleles;$i++){
				$multi_alleles[$i]=~s/[\r]//;
				my @annot_info=split("\t",$annot);my @addl_orig_info=split("\t",$orig);
				my $addon='';
				if ($multi_alleles[$i] eq '.'){
					next ALLELE;
				}elsif ($annot_info[$idx{'ANVR_TYPE'}]!~/(exonic|spli)/i && $i>0){
#					next ALLELE if $i>0;
				}
				if ($opt_f=~/vcf/i){
					# for (my $i=0;$i<=$#annot_info;$i++){
					# 	print "$i\t$annot_info[$i]\n";
					# }
					# <STDIN>;
					my $fmtidx=0;
					if($annot=~/^(.*$$aref[8]:SMAF)/){
						my $foundspot=$1;$$aref[8].=":SMAF" if ($$aref[8]!~/SMAF/);
						$fmtidx++ while ($foundspot=~/\t/g);
					}elsif($annot=~/SMAF/ && $annot=~/^(.*$$aref[8])/){
						my $foundspot=$1;
						$fmtidx++ while ($foundspot=~/\t/g);
					}
					chop $$aref[7] until ($$aref[7]!~/[;,:]$/ || $$aref[7] eq '');
					$$aref[7].=";" if ($$aref[7]!~/[;,]$/);
					if ($i==0 || $opt_w){
						# print Dumper (\%idx);die;
						# print Dumper (\%fdi_headers);die;
						foreach my $id (sort keys %idx){
							
							next if ($id=~/^[a-z\_]*$/);
							next if ($id=~/(Query (Start|end)|Allele\s*\d|^chr$)/i);
							# if ($debug) {print "$orig_info[1]:working on $multi_alleles[$i] and $id($i)...\n$last_annot\n";<STDIN>;}
							# $addon= ($id=~/ANNOT/ && $#multi_alleles>0 && $annot_info[$idx{'ANVR_TYPE'}]=~/(exon|splic)/)?"($multi_alleles[$i])":'' ;#changed this to next line
							$addon='' if ($i==0);
							# print "[DEBUG] working on $id\n";
							# print "what's this $id($i)?($addon)\n";
							$addon= ($id!~/(ANVR_GENE|ANVR_TYPE)/i)?"($multi_alleles[$i])":'' ;
							# if ($debug){print "\t$addon for $id...\n";}
							my $val=$annot_info[$idx{$id}];
							$val=~s/\s/_/g;
							$val=~s/;/:/g if ($id!~/;/i);
							$val=~s/Score=([\d\.]+).Name=(\S+)/\2\(\1\)/g;#$val=~s/Score=(\d\.+)/\1/;
							$val=~s/Score=([\d\.]+)/\1/g;#$val=~s/Score=(\d\.+)/\1/;
							$val=~s/Name=//g;
							# print "working on $id\n";<STDIN>;
							my $tmpannot='';
							if ($FunSeq2 && $id=~/;/ ){
								my @FunSeq2_arr=split(";",$val);
								for(my $j=0;$j<=$#FunSeq2_headers;$j++){
									# $last_annot.=uc("FunSeq2_".$FunSeq2_headers[$j])."$addon=$FunSeq2_arr[$j];" and "adding\n" if ($FunSeq2_arr[$j] ne '.');##20160104
									$last_annot=addAnnot(uc("FunSeq2_".$FunSeq2_headers[$j]), $FunSeq2_arr[$j], $addon,$last_annot);
								}
							}elsif($exac && $id=~/ExACv/i && $val=~/:/){
								my @exac_arr=split(":",$val);
								foreach my $foo (keys %exac){
									# $last_annot.= uc$foo."$addon=$exac_arr[$exac{$foo}[1]];";##20160104
									$last_annot=addAnnot(uc$foo,$exac_arr[$exac{$foo}[1]],$addon,$last_annot);
								}
							}elsif($val=~/=\S+/ && $id=~/(exac(\d+)|EXOME AGGREGATION.*?(V\d+))/i ){
								my $ver=($3)?$3:"v$2";
								$val=~s/=/$ver$addon=/;
								$last_annot.="$val;";
							}elsif (exists $fdi_headers{$id}){
								# $last_annot.="$fdi_headers{$id}$addon=$val;" if ($annot_info[$idx{$id}]!~/(^\-|ERR|^)$/i);#/mnt/webrepo/fr-s-abcc-avia1/public/upload/20140804170912.IG.AviaUpload.vcf
								$last_annot=addAnnot($fdi_headers{$id},$val,$addon,$last_annot);
							}else{
								# $last_annot.="$id$addon=$val;" if ($annot_info[$idx{$id}]!~/(^\-|ERR|^)$/i);#ANVR_GENE=$exinfo[1];ANVR_ANNOT=$exinfo[2];";
								# print "about to add $id and $val to $last_annot\n\n";
								$last_annot=addAnnot($id,$val,$addon,$last_annot);
								# print "Outside of addAnnot\n$last_annot\n\n";
							}
						} 
						# print $last_annot; <STDIN>;
					}elsif($i>0 && $annot_info[$idx{'ANVR_TYPE'}]=~/(exon|splic)/){
						my $val=$annot_info[$idx{ANVR_ANNOT}];$val=~s/\s/_/g;
						$last_annot.="ANVR_ANNOT($multi_alleles[$i])=$val;";
						if ($debug){print "about to print annovar annotation to the last annot!\n";}
					}
					if ($i!=$#multi_alleles){
						# $debug=1;print "\t[DEBUG]advancing annot,$#multi_alleles,$i\n\t".join ("\t",@{$aref})."\n";<STDIN>;
						$annot=<ANNOT>;
						$annot=~s/[\n\r]//g;
					}
					if ($fmtidx){
						for (my $d=9;$d<=$#{$aref};$d++){#these are the sample rows
							$$aref[$d]=$annot_info[$fmtidx+$d-8];
						}
					}
					if ($debug){print "Adding this to aref[7]...$last_annot";<STDIN>;}
					
				}else{
					my @arr=($orig_info[0],$orig_info[1], getIDField($ct,''),$orig_info[3],$orig_info[4],'.','.');
					foreach my $id (sort keys %idx){
						next if ($id=~/^[a-z\_]*\_\([ATGC]*\){0,1}$/);#ignore anything with lower case as it is just info
						$arr[7].="...$id$addon=$annot_info[$idx{$id}];" if ($annot_info[$idx{$id}]!~/(^\-|ERR)$/i);#ANVR_GENE=$exinfo[1];ANVR_ANNOT=$exinfo[2];";
					} 
					push (@arr,"GT:GC");
					
					if ($#orig_info>14){
						push (@arr,@orig_info[5..$#orig_info]);
					}
		#			push (@orig_info,join ("\t",@arr,"\n"));
					
				}
			}
			$$aref[7].=$last_annot  ;
			chomp $$aref[$#{$aref}];

			# if ($debug){print "LN".__LINE__.":$#multi_alleles|$annot\n======\n". join ("\t",@{$aref})."\n";<STDIN>;}
			chop $$aref[7] until ($$aref[7]!~/;$/ || $$aref[7] eq '');
			$$aref[7]=~s/;{2,}/;/g;$$aref[7]=~s/^[\.\s]*;\s*//;$$aref[7]=~s/[\.,:];/;/g;$$aref[7]=~s/[,;:]$//;
			print FILE join ("\t",@{$aref})."\n";
		}else{
			my @annot_info=split("\t",$annot);
			my @arr=($orig_info[0],$orig_info[1], getIDField($ct,''),$orig_info[3],$orig_info[4],'.','.');
			foreach my $id (sort keys %idx){
				next if ($id=~/^[a-z\_]*\_\([ATGC]*\){0,1}$/);#ignore anything with lower case as it is just info
				my $val=$annot_info[$idx{$id}];$val=~s/\s{1,}/_/g;
				$arr[7].="$id=$val" if ($annot_info[$idx{$id}]!~/(^\-|ERR)$/i);#ANVR_GENE=$exinfo[1];ANVR_ANNOT=$exinfo[2];";
			} 
			push (@arr,"GT:GC");
			
			if ($#orig_info>14){
				push (@arr,@orig_info[5..$#orig_info]);
			}
			print FILE join ("\t",@arr)."\n";
#			push (@orig_info,join ("\t",@arr,"\n"));
		}
		
	}
  
  }
#}else{
#	#just put headers and delete one column
#	my ($header,$footer)=getAuxHeaders();
#	#printing VCF header...
#	print FILE "$header". findDBheaders()."$footer";#this does not have a new line at the end
#	#need to find out if there is patient information in the file...ie. columns 7-n in the file
#	#txt input files do not have any headers
#	die "Your file does not exist or is empty\n" if (!-e "$annot_fn" || -z $annot_fn);
#	my  @info;#holder for batched writing;
#	open (P,"<$opt_i") or die "Cannot open $opt_p for reading\n";
#	my $ct=0;
#	while (my $line=<P>){
#		chomp $line;$ct++;
#		
#		$line=~s/[\r\n]//g;
#		if ($ct%$FILE_LOAD==1  ){
#			print FILE join ("",@info);
#			@info=();#clear info;
#			doStuff($ct);
#		}
#		next if $ct==1;#header
#		my @inparr=split("\t",$line);#user input
#		my @arr=($inparr[4],$inparr[5], getIDField($ct,''),$inparr[7],$inparr[8],'.','.',"ANVR_TYPE=$inparr[0];ANVR_GENE=$inparr[1];ANVR_ANNOT=$inparr[2];");
#		$arr[7].="$hash{$ct+1}" if (exists $hash{$ct+1});
#		if ($inparr[9] ne ''){push (@arr,"COMMENT=$inparr[9]");}
#		push (@arr,"\tGT:GC");
#		
#		if ($#inparr>9){
#			push (@arr,@inparr[10..$#inparr]);
#		}
#		push (@info,join ("\t",@arr,"\n"));
#	}
#	close P;
#	if ($#info>-1){
#		print FILE join ("",@info);
#	}
#}

close FILE;
	
	
################################# SUBROUTINES #######################################
sub addAnnot(){
	my $xheader=shift;
	my $xannot=shift;
	my $xaddon=shift;
	my $curr=shift;
	my $xannotx;
	if (!$opt_w){
		$xaddon='';
	}
	if ($curr!~/$xheader=[^=]*$xannot/ && $xannot!~/(^\.$|ERR|^\-$);*/ ){
		# if ($debug){print "\tadding\n\t\t $xannot for $xheader ($xaddon)\n";}
		$curr.="$xheader=$xaddon$xannot;";
	}elsif ($xaddon && $xannot!~/(^\.$|ERR|^\-$)/ && $curr=~/$xheader=.*$xannot/){
		if ($curr=~/$xheader=\(([ATCG\-,]+)\)$xannot/){
			my $other=$1;
			# print "start with $xaddon\n";
			$xaddon=~s/(\w+)/$other,$1/;
			# print "$curr\n\n replace $xaddon??($xheader)\n";
			$curr=~s/$xheader=\($other\)$xannot/$xheader=$xaddon=$xannot/;
			# print "$curr\n";<STDIN>;
		}else{
			print "[DEBUG] $curr not have any allele!\n";
		}
	# }elsif ($curr=~/$xheader.=$xannot/ || $xannot=~/(^\.$|ERR|^\-$)/ || !$xaddon){
	# 	print "[DEBUG] Skip because exists or null ...$xannot\n";
	# }else{
	# 	print "[DEBUG] skipping anything for \n\tHEADER:$xheader($xaddon)\n\tXANNOT:$xannot\n\tCURR:$curr\n";<STDIN>;
	}
	return $curr;
}
sub getIDField{
	my $id=shift;my $curr=shift;
	if ($opt_t=~/(rsid|append|none)/i || ($opt_t=~/user_defined/ && $opt_n)){
		my $replaced_name;
		if ($debug){print "LN".__LINE__.":looking for $id...$$name_href{$id+$OffbyOne}\n".Dumper (\%{$name_href});<STDIN>;}
		switch ($opt_t){
			case 'rsid'	{$replaced_name="MOO2".$$name_href{$id+$OffbyOne};}
			case 'append'	{$replaced_name=($curr ne '.')?$curr.";$$name_href{$id+$OffbyOne}" : $$name_href{$id+$OffbyOne};}
			case 'none' {$replaced_name="ABCC_".sprintf("%04d",$id);}
			case 'user_defined'	{$replaced_name="MOO1".$$name_href{$id+$OffbyOne};}
		}
		chop $replaced_name if ($replaced_name=~/;$/);
		$replaced_name= ($replaced_name=~/^\s{0,}$/)?'.':$replaced_name;
		$curr=$replaced_name;
	}elsif ($curr eq ''){
		$curr='.';
	}
	return "$curr";
}
sub findDBheaders{
	my $text='';
	$text.="##INFO=<ID=ANVR_GENE,Number=.,Type=String,Description=\"Gene Annotations from ANNOVAR\">\n".
		"##INFO=<ID=GENE_ORI,Number=.,Type=String,Description=\"Gene Orientation (if exists) from ANNOVAR\">\n".
		"##INFO=<ID=ANVR_TYPE,Number=.,Type=String,Description=\"Genic Type Annotation from ANNOVAR\">\n".
		"##INFO=<ID=ANVR_ANNOT,Number=.,Type=String,Description=\"Annovar Annotations. All annotations correspond to values for first variant in 'ALT' Col (multi-variants ignored)\">\n";
	if (scalar(@db_searched_files)){
		foreach my $searchfile (@db_searched_files){
			next if (defined $opt_n && $searchfile=~/$opt_n/i);
			my $file=`basename $searchfile`;chomp $file;
			$file=uc($file);
			$file=~s/($annovar_base\.|hg1[89]\_|mm[8|9]\_|.txt)//gi;
			my $tmpinfo;
			if (-e "avia.database.info.txt"){
				if ($file=~/FunSeq$/i){
					my %fs_hash=("FunSeq2_altcds"=>"FunSeq2 annot for coding/noncoding (yes/no)",
						"FunSeq2_variant.annotation.cds"=>"FunSeq2 annot for coding variants",
						"FunSeq2_network.hub"=>"FunSeq2 score for whether the gene is a hub",
						"FunSeq2_gene.under.negative.selection"=>"FunSeq2 boolean for whether gene is under negative selection",
						"FunSeq2_ENCODE.annotated"=>"FunSeq2's ENCODE annotations",
						"FunSeq2_motif.breaking"=>"FunSeq2 flag for whether it breaks a motif for binding",
						"FunSeq2_sensitive"=>"FunSeq2 score for sensitivity",
						"FunSeq2_ultra.sensitive"=>"FunSeq2 score for ultrasensitivity",
						"FunSeq2_target.gene"=>"FunSeq2 annotation for Gene(s)",
						"FunSeq2_coding.score"=>"FunSeq2 Prioritization Score (if coding)",
						"FunSeq2_noncoding.score"=>"FunSeq2 Prioritization Score (if non-coding)");
					foreach my $key (keys %fs_hash){
						$text.="##INFO=<ID=$key,Number=.,Type=String,Description=\"Annot from AVIAv2.0:$fs_hash{$key}\">\n";
						# print "##INFO=<ID=$key,Number=.,Type=String,Description=\"Annot from AVIAv2.0:$fs_hash{$key}\">\n";<STDIN>;
					}
				}elsif ($file=~/FunSeq2/i && -e "FunSeq2.headers.out"){
					open (FS,"<FunSeq2.headers.out");
					my @fs_arr=split(/;/,<FS>);
					print "Using FunSeq2.headers.out\n\n";
					close FS;chomp $fs_arr[$#fs_arr];
					foreach my $key (@fs_arr){
						$key=uc("FunSeq2_".$key);
						$text.="##INFO=<ID=$key,Number=.,Type=String,Description=\"Annot from AVIAv2.0:FunSeq2 Annotations\">\n";
						# print "working on $key and $fs_arr[$key]\n";
						# $fdi_headers{$key}=uc($fs_arr[$key]);
						$FunSeq2++;
					}
					# print Dumper (\%fdi_headers);die;
				}elsif($file=~/ExACv\d+(\_\d+)*/i){
					my $ver=$1;
					$ver||="";
					%exac=("ExACv3_AltAC_Hom"=>["ExACv3$ver Alt Count",0],
						"ExACv3_All_MAF"=>["ExACv3$ver All Alleles MAF" ,1],
						"ExACv3_NFE_MAF"=>["ExACv3$ver Non-Finnish MAF",2],
						"ExACv3_Highest_MAF"=>["Highest MAF",3],
						"ExACv3_Ethnicity"=>["ExACv3$ver Highest MAF Ethnicity",4]);
					foreach my $key (keys %exac){
						$text.="##INFO=<ID=".uc($key).",Number=.,Type=String,Description=\"Annot from AVIAv2.0:$exac{$key}[0]\">\n";
						$fdi_headers{uc($key)}=uc($exac{$key}[0]);
						# print "##INFO=<ID=$key,Number=.,Type=String,Description=\"Annot from AVIAv2.0:$exac{$key}\">\n";<STDIN>;
					}
					$exac++;
				}elsif ($tmpinfo){
					$tmpinfo=`grep -i $file "avia.database.info.txt" -m 1 `;chomp $tmpinfo;
					# print "found $tmpinfo\n";
					$tmpinfo=~s/[\r]//g;
					$tmpinfo=~s/&lt;/</g;
					$tmpinfo=~s/&gt;/>/g;
					$tmpinfo=~s/&quot;/'/g;
					my @arry=split("\t",$tmpinfo);
					$fdi_headers{$arry[1]}=uc($arry[0]);
					$fdi_headers{$arry[0]}=uc($arry[0]);
					$tmpinfo=":Description=\"$arry[1]\":Source=$arry[2]:Version=$arry[3]";
				}
				$text.="##INFO=<ID=$file,Number=.,Type=String,Description=\"Annot from AVIAv2.0$tmpinfo\">\n";
				
			}else{
				$text.="##INFO=<ID=$file,Number=.,Type=String,Description=\"Annot from AVIAv2.0\">\n";
			}
		}
		# print Dumper (\%fdi_headers);die;
	}
	return $text;
}

sub getAuxHeaders{
	my $date=`date "+%Y%m%d"`;chomp $date;
	my $h=qq(##fileformat=VCFv4.0
##fileDate=$date
##source=user_defined/ABCC Impact Analysis Pipeline
);
	my $f=join("\t","#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT");
	return ($h,$f);
}
sub readAnnotBlockFDI{
	my $fn=shift;
	my $startln=shift;
	my $endln=shift;
	my $arr=shift;
	my %tmphash;
	open (FH,"<$fn") or die "cannot open $fn at LINE".__LINE__."\n";
	my $ct=0;my %idx;my @arr;
	if (ref($arr)=~/ARRAY/){
	}else{
		@arr=split("\t",`head -n1 $fn`);chomp $#arr;
		$arr=\@arr;
	}
	for (my $i=0;$i<=$#{$arr};$i++){
		if ($$arr[$i]=~/^(#){0,1}ANNOVAR/i){
			$idx{'ANVR_TYPE'}=$i;
		}elsif ($$arr[$i]=~/#(\S+)/i){
			my $foo=$1;
			if ($$arr[$i]=~/;/){
				$idx{"FUNSEQ2"}=$i;
				$idx{uc($foo)}=$idx{"FUNSEQ2"};
			}else{
				$idx{uc($foo)}=$i;
			}
		}elsif ($$arr[$i]=~/^Gene$/i){
			$idx{'ANVR_GENE'}=$i;
		}elsif ($$arr[$i]=~/(Rel pos|Annot Feat)/i){
			$idx{'ANVR_ANNOT'}=$i;
		}
	}
	while (<FH>){
		my @exinfo=split("\t",$_);
		$ct++;
		if ($ct>=$startln && $ct<=$endln+1){
			if (exists $tmphash{$ct}){die "this shouldn't already exist!\n";
				$tmphash{$ct}.=",$exinfo[2]";
			}else{
				foreach my $id (sort keys %idx){
					$tmphash{$ct}.="$id=$exinfo[$idx{$id}];" if ($exinfo[$idx{$id}] ne '-');#ANVR_GENE=$exinfo[1];ANVR_ANNOT=$exinfo[2];";
				}
			}
		}
		last if ($ct>$endln);
	}
	close FH;
	return (\%tmphash,$arr);
}
