#!/usr/local/bin/perl 
#use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
umask(0000);
use Data::Dumper;
=head
This module describes ABCC/LMT specific modules for additional calculations not done by ANNOVAR
v2.0
bug fixes to the output where errors occur
=cut
sub SAFE_PROCESS{
	my $cmd=shift;
	my $err_msg=shift;
	my $die=shift;
	print "$cmd\n" if ($die);
	eval {
		system ("$cmd\n");
	};
	if ($@){
		print "$err_msg($@)\n";
		die $err_msg if ($die);
	}else{
		return 1;
	}
}
sub _index{
	my $db_dir=shift;
	my $db_file=shift;
	my $name_col=shift; #optional
	return if (!$db_dir || !$db_file);
		##assume that the format is "chr\tstart\tstop\tcomments\n";
	if (SAFE_PROCESS ("cp $db_dir/$db_file $db_dir/$db_file.TMP\n","Could not copy at line ".__LINE__."\n",0)){
		#try to figure out the format of the file
		my @head;my %types=();
		my $orig_size=-s "$db_dir/$db_file";#keep track of original size in case we need it;
		@head=`grep -ve '^#' $db_dir/$db_file.TMP | head -n 10` if (-e "$db_dir/$db_file.TMP");
		if ($#head>0){
			my $filetype;
			foreach my $line (@head){
				my @info=split("\t",$line);
				my $regex='(chr){0,1}[\dXYMT]{1,2}\t\d+\t\d+\t';
				if ($line=~/^$regex/){
					$filetype='A';
				}elsif ($line=~/^\d+\t$regex/ && $info[1]=~/(chr){0,1}[\dXYMT]+/){
					$filetype='B';
				}elsif ($line=~/^\w+\t\w+\t$regex/ && $info[2]=~/(chr){0,1}[\dXYMT]+/){
					$filetype='C';
				}else{
					print "cannot determine filetype from $regex\n\t$line\n";
				}
				last if ($filetype ne '');
			}
			if (!SAFE_PROCESS("perl $bin/index_annovar_ABCC.pl $db_dir/$db_file.TMP --outfile $db_dir/$db_file --filetype $filetype","indexing failed",0)){
				_testIndexing("$db_dir/$db_file",$orig_size);die "for testing";
				return 0;
			}
			
		}else{
			_testIndexing("$db_dir/$db_file",$orig_size);
			warn "could not find the header in this file $db_dir/$db_file\n";
			return 0;
			#do not index and leave this routine
		}
	}else{die;
		return 0;
	}
	$name_col||=3;#default is 3
	my $basename=$db_file;
	$basename=~s/^(hg1\d|mm\d+)\_([\w\-\.]+).txt/$2/g;
	if ( -e "$db_dir/lkup.txt"){
		my $test_exists=`grep -P '^$basename' $db_dir/lkup.txt`;chomp $test_exists;
		if (!$test_exists || $test_exists!~/$basename\t\S*\t[\d\,]+\t\S*\t$name_col/){
			open (LKUP,">>$db_dir/lkup.txt") or die "cannot open lkup.txt in $db_dir for appending\n$!";
		}else{
			return 1;
		}
	}else{
		open (LKUP,">$db_dir/lkup.txt") or die "cannot open lkup.txt in $db_dir for writing\n$!";
	}
	print LKUP "$basename\t\t0,1,2\t\t$name_col\n";
	close LKUP;
	unlink ("$db_dir/$db_file.TMP");
	return 1;
}
sub _testIndexing{
	my $xbase=shift;my $xsize=shift;
	if (-s "$xbase.TMP" == $xsize){
		system ("mv $xbase.TMP $xbase\n");
		warn "could not index $xbase ($@)($!)\n";
	}elsif (-s "$xbase" == $xsize){
		system ("rm $xbase.TMP\n") if (-e "$xbase.TMP");
	}else{
		die "Error processing the indexing phase!!!\n";
	}
}
sub _getABCC_DBInfo{
	my $db_loc=shift;
	my $xdbname=shift;
	print "in _getABCC_DBInfo...\n";
	if ( -e "$db_loc/lkup.txt"){
		my $info=`grep -P '^$xdbname\t' $db_loc/lkup.txt`;chomp $info;
		if ($info eq ''){
			print "Found ($info)\n";
			@info=($xdbname,'','0,1,2','','3');
		}else{
			print "FOUND:$info!!!\n";
			@info=split("\t",$info);
		}
	}else{
		@info=($xdbname,'','0,1,2','','3');
		return (\@info);
	}
	return (\@info);
}
sub _calculateAdditionInfo{
	my $trx_len_href=shift;
	my $info=shift;
	my ($name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend, $name2) = @$info;
	my @exonstart=@$exonstart;
	my @exonend=@$exonend;
	$gene_ori{$name2}=$dbstrand if (! exists $gene_ori{$name2});
	my ($cdsoffset,$trx_lastCDS,$lastCDS_enum,$isE1UTR,$isLastUTR);
	if ($dbstrand eq '-'){#validated with WDR1 MD1081 rs13441
		#query        				 ----
		#gene <------>   <------UUUUUU>  <UUUUUU>  U=UTR
		my $last=0;
		for (my $i=$#exonstart;$i>=0;$i--){
			#lastUTR_len is the length of the TOTAL full UTRs on the 3' end.
			if ($cdsstart>=$exonstart[$i] && $cdsstart<=$exonend[$i]){
				$lastCDS_enum=$i;$trx_lastCDS+=(abs($exonstart[$i]-$exonend[$i])+1);
				last;
			}else{
				$trx_lastCDS+=(abs($exonstart[$i]-$exonend[$i])+1);
			}
		}
		if ($cdsstart==$cdsend || ($cdsstart==$exonstart[0] && $cdsend==$exonend[$#exonend])){
			$cdsoffset=0;
			$lastCDS_enum=0;
		}else{
			$isE1UTR=( ($cdsend>$exonstart[$#exonstart]) && ($cdsend<$exonend[$#exonstart]))?1:0;
			$isLastUTR=( $cdsstart>$exonstart[0] && $cdsstart<$exonend[0])?1:0;#cds bound is within the first exon (i.e. last exon)
			if ($isE1UTR){
				$cdsoffset=abs($cdsend-$txend);#?do i need to add 1?
			}else{#validated KRAS MD1055 rs4362222  
				$cdsoffset=0;my $index=$#exonstart;
				until (($cdsend>=$exonstart[$index] && $cdsend<=$exonend[$index]) || $index==0){
					$cdsoffset+=abs($exonstart[$index]-$exonend[$index])+1;
					$index--;
				}
				$cdsoffset+=abs($exonend[$index]-$cdsend);
			}
		}
	}elsif ($dbstrand eq '+'){#tested CYBB MD1012 rs2228117 
		#query  				----
		#gene <UUUUU>   <UUUUUU----> <--->
		for (my $i=0;$i<=$#exonstart;$i++){
			if ($exonstart[$i]<=$cdsend && $exonend[$i]>=$cdsend){
				$lastCDS_enum=$i;$trx_lastCDS+=(abs($exonstart[$i]-$exonend[$i])+1);
				last;
			}else{
				$trx_lastCDS+=(abs($exonstart[$i]-$exonend[$i])+1);
			}
		}
		if ($cdsstart==$cdsend || ($cdsstart==$exonstart[0] && $cdsend==$exonend[$#exonend])){
			$cdsoffset=0;
			$lastCDS_enum=0;
		}else{
			$isE1UTR=($cdsstart>$exonstart[0] && $cdsstart<$exonend[0])?1:0;
			$isLastUTR=( $cdsend>$exonstart[$#exonstart] && $cdsend<$exonend[$#exonend])?1:0;
			if ($isE1UTR){#validated MD998 G6PC3 rs3815076 
				$cdsoffset=abs($cdsstart-$txstart);#do I need to add 1?
			}else{#tested with rs2275774 C10orf18 no MD  
				$cdsoffset=0;my $index=0;
				until (($cdsstart>=$exonstart[$index] && $cdsstart<=$exonend[$index]) || $index>=$#exonstart){
					$cdsoffset+=abs($exonstart[$index]-$exonend[$index])+1;
					$index++;
				}
				$cdsoffset+=abs($exonstart[$index]-$cdsstart);
			}
		}
	}else{
		die "cannot find strand at LN".__LINE__."\n";
	}
	die "($trx_lastCDS) cannot be smaller that last utr len $lastCDS_enum for \t".join(",",@$info,"\n",@exonstart,"\n",@exonend)."\n" if ($lastCDS_enum>$trx_lastCDS && $name!~/^NR/i);
	$$trx_len_href{$name}="$cdsoffset;$trx_lastCDS;$lastCDS_enum;$isE1UTR;$isLastUTR";
	return ($trx_len_href);
}
sub _newprocessNextQueryBatchByGene_withExtras{
	my $debug=0;
	my ($queryfh, $batchsize, $genedb, $geneidmap, $cdslen, $mrnalen) = @_;
	my (%refseqvar);
	my ($chr, $start, $end, $ref, $obs);
	my ($name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend, $name2, $cdsarr);
	my ($invalid);
	my ($linecount, $invalidcount) = qw/0 0/;
	for my $i (1 .. $batchsize) {					#process up to batchsize variants
		my $nextline = <$queryfh>;				#read the next line in variant file
		defined $nextline or last;
		$nextline =~ s/[\r\n]+$//;
		$linecount++;						#linecount does not include the comment line		
		if ($nextline =~ m/^#/ and $comment) {			#comment line start with #, do not include this is $linecount
			print OUT "$nextline\n" if (!$noheader);
			next;
		}elsif ($nextline=~m/^#/ ){
			next;
		}
		
		$invalid = 0;
		
		my @nextline = split (/\s+/, $nextline);
		($chr, $start, $end, $ref, $obs) = @nextline[@avcolumn];
		if ( not (defined $chr and defined $start and defined $end and defined $ref and defined $obs ) and !$novar) {
			$invalid++;
		} else {
			($ref, $obs) = (uc $ref, uc $obs);
			$zerostart and $start++;
			$chr =~ s/^chr//;
			if ($chr =~ m/[^\w]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {
				$invalid++;
			}elsif ($novar){#added if you just wanted to see where the region fell in the gene annotations (ie. annotating nonB)#hv added 20120329
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
			print INVALID $nextline, "\n";			#invalid record found
			print OUT "ERR\tERR\tERR\tERR\t$nextline\n" if ($keepline);
			$invalidcount++;
			next;
		}elsif($keepline){
			print NEW "$nextline\n";
		}
# print __LINE__. ":valid $nextline\n";<STDIN>;
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
##hv added block
		my (%seen,%trx_len_hash,$trx_len_href);  #hv added this!!
		$trx_len_href=\%trx_len_hash;
##end block
		for my $nextbin ($bin1 .. $bin2) {
			%trx_len_hash=();#clear hash from last set #hv added
			exists $genedb->{$chr, $nextbin} or next;		#this genome bin has no annotated gene (a complete intergenic region)
			for my $nextgene (@{$genedb->{$chr, $nextbin}}) {	#when $genedb->{$chr, $nextbin} is undefined, this automatically create an array!!!
				($name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend, $name2) = @$nextgene;
				#interesting examples: snippets from hg19_refGene.txt
				#624     NM_001005160    chr11   -       5152921 5153872 5152921 5153872 1       5152921,        5153872,        0       OR52A5  cmpl    cmpl    0,
				defined $name2 or print STDERR "WARNING: name2 field is not provided for transcript $name (start=$txstart end=$txend)\n" and $name2='';
				$seen{$name, $txstart} and next;		#name and txstart uniquely identify a transcript and chromosome position (sometimes same transcript may map to two nearby positions, such as nearby segmental duplications)
				$seen{$name, $txstart}++;			#a transcript may be in two adjacent bins, so if one is already scanned, there is no need to work on it again
				my $current_ncRNA;
#hv added next block
				print "working on $name,$name2 at LN".__LINE__."\n" if ($debug);
				#first determine if the 5' UTR is in the first exon
				#then calculate the cdsoffset
				my $cdsoffset= my $isE1UTR=my $isLastUTR= my $lastCDS_enum=my $trx_lastCDS=0;#trx_lastCDS is a variable containing the mrna pos up to the last coding exon.  We do this so we can better calculate the relative position into the 3UTR
				if (!exists $$trx_len_href{$name}){
					($trx_len_href)=_calculateAdditionInfo($trx_len_href,$nextgene);
				}
				($cdsoffset,$trx_lastCDS,$lastCDS_enum,$isE1UTR,$isLastUTR)=split(";",$$trx_len_href{$name});
				##############################ends hv block
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
					my $distToNext=abs($txstart-$end)-1;
					if ($end > $txstart - $neargene) {
						if ($dbstrand eq '+') {
							$upstream{$name2}="-$distToNext;" if (!exists $upstream{$name2});
						} else {
							$downstream{$name2}="+$distToNext;" if (!exists $downstream{$name2});
						}
					} else {
						last;						#if transcript is too far away from end, end the search of the bins
					}
				} elsif ($start > $txend) {
					#query            ---
					#gene  <-*----*->
					if (not $foundgenic and $start < $txend + $neargene) {
						if ($dbstrand eq '+') {
							$downstream{$name2}="+".(abs($txend-$start)-1).";" if (!exists $downstream{$name2});
						} else {
							$upstream{$name2}="-".(abs($txend-$start)-1).";" if (!exists $upstream{$name2});
						}
					}
#				} elsif ($cdsstart == $cdsend+1) {				#non-coding RNA (could be microRNA, or could be due to lack of CDS annotation for mRNA such as NR_026730 or BC039000). Previously we already did cdsstart++ so here the cdsstart is more than cdsend
#					if ($start >= $txstart and $start <= $txend or $end >= $txstart and $end <= $txend or $start <= $txstart and $end >= $txend) {
##						$ncrna{$name2}++;
#						$ncrna{$name2}.="$name;";
#						$foundgenic++;
#					}
				} else {							#query overlaps with coding region of gene
					if ($cdsstart == $cdsend+1) {				
						#non-coding RNA (could be microRNA, or could be due to lack of CDS annotation for mRNA such as NR_026730 or BC039000). Previously we already did cdsstart++ so here the cdsstart is more than cdsend
						if ($start >= $txstart and $start <= $txend or $end >= $txstart and $end <= $txend or $start <= $txstart and $end >= $txend) {
							$ncrna{$name2}.="$name;";
							$foundgenic++;
						}
						
						#now treat this ncRNA as if it is a protein-coding gene
						($cdsstart, $cdsend) = ($txstart, $txend);
						$current_ncRNA++;		#current transcript is a noncoding transcript
					}
					my ($lenintron, $lenexon) = (0, 0);			#cumulative intron/exon length at a given exon
					my ($rcdsstart, $rvarstart, $rvarend);			#start of coding and variant in reference mRNA sequence
					my $foundexonic;
					my @exonstart = @$exonstart;
					my @exonend = @$exonend;
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
								$splicing{$name2}.="$name," if ($splicing{$name2}!~/$name,/);		#when query start site is close to exon start or exon end
							}
							if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1 or $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
								$splicing{$name2}.="$name," if ($splicing{$name2}!~/$name,/);		#when query end site is close to exon start or exon end
							}
							if ($start <= $exonstart[$k] and $end>=$exonstart[$k] or $start <= $exonend[$k] and $end >= $exonend[$k]) {
								$splicing{$name2}.="$name," if ($splicing{$name2}!~/$name,/);		#when query encompass the exon/intron boundary
							}
							
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if this splicing leads to amino acid change (rather than UTR change)
							if ($splicing{$name2} and $start==$end and $start>=$cdsstart) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$lenexon -= ($exonend[$k]-$exonstart[$k]);		#formerly "$exonend[$k]-$exonstart[$k]+1"; changed this line 2011oct01 given German's comments. The coding portion for acceptor site should be added by one, but for donor site should be fine
									$splicing_anno{$name2} .= "$name:c.$lenexon-" . ($exonstart[$k]-$start) . "$ref>$obs,";
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
									my ($diff,$mpos);
									#here the trick begins to differentiate UTR versus coding exonic
									if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
										#query  ----
										#gene     <--*---*->
										#original code	$utr5{$name2}++;		#positive strand for UTR5\
										#hv note: no UTR5 in GRID_Human_37.I_SnpImpacts so cannot really test
										$diff=abs($end-$cdsstart);
										$mpos=$cdsoffset-$diff;
										$utr5{$name2}.="$name:E".sprintf("%02d",$k+1).":m.$mpos:c.-$diff;";
									} elsif ($start > $cdsend) {#positive strand for UTR3
										#query             ----
										#gene     <--*---*->
										#hv added next block:
										$diff=0;$mpos=0;
										if ($end>$exonstart[$lastCDS_enum] && $end<$exonend[$lastCDS_enum]){
											#it is in the last coding exon so just calculate the distance from the cds pos
											$diff=(abs($cdsend - $end)+1);$mpos+=abs($cdsend-$exonstart[$lastCDS_enum])+1;
										}else{
											#need to find the last coding exon UTR_length
											$diff=abs($exonend[$lastCDS_enum]-$cdsend);#calculates the UTR of the last cds/utr exon;
											$mpos+=abs($cdsend-$exonstart[$lastCDS_enum]);
											for (my $z=$lastCDS_enum+1;$z<$k;$z++){
												$diff+=abs($exonstart[$z]-$exonend[$z]);
											}
											#add the current distance
											$diff+=abs($end-$exonstart[$k])+1;
										}
										print "\t$diff and $trx_lastCDS\n";
										$mpos+=$trx_lastCDS+ $diff;
										#end hvblock
										$utr3{$name2}.="$name:E".sprintf("%02d",$k+1).":m.$mpos:c.+$diff;";
#										$utr3{$name2}.="$name:E".sprintf("%02d",$k+1).":+".abs($start-$cdsend).";";#hvmodh
#										$utr3{$name2}++;		
									} else {							
										$exonic{$name2}.="$name:o$cdsoffset;" ;
										not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '+', $i, $k+1, $nextline];	#refseq CDS start, refseq variant start. obs is non-zero (obs is specified by user)
									}
									$foundgenic++;
									last;
								} elsif ($k and $start > $exonend[$k-1]) {	#intronic on the positive
									#HUEFLAG	
									#calculate if it is closer to the exon$k or exon$k+1
									print "at LN".__LINE__." looking at intron...\n" if ($debug);
									my ($flg,$distToExon)=_getNearest($start,$exonend[$k-1],$exonstart[$k],'+');
									print "LN".__LINE__. ":got $k:$flg,$distToExon and $k+$flg \n" if ($debug);
									$intronic{$name2}.="$name:E".sprintf("%02d",($k-$flg+1)).":$distToExon;";#End hueflag 
									#original code: 	$intronic{$name2}++
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
								my ($diff,$mpos);
								#here is the trick begins to differentiate UTR versus coding exonic
								if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
									#query  ----
									#gene     <--*---*->
									$diff=abs($end-$cdsstart);#hv
									$mpos=$cdsoffset-$diff;#hv
									$utr5{$name2}.="$name:E".sprintf("%02d",$k+1).":m.$mpos:c.-$diff;";#hv
#									$utr5{$name2}++;		#positive strand for UTR5 #original
								} elsif ($start > $cdsend) {
									#query             ----
									#gene     <--*---*->
									#hv block starts
									$diff=0;$mpos=0;
									if ($end>$exonstart[$lastCDS_enum] && $end<$exonend[$lastCDS_enum]){
										#it is in the last coding exon so just calculate the distance from the cds pos
										$diff=(abs($cdsend - $end));$mpos+=abs($cdsend-$exonstart[$lastCDS_enum]);
									}else{
										#need to find the last coding exon UTR_length
										$diff=abs($exonend[$lastCDS_enum]-$cdsend)+1;#calculates the UTR of the last cds/utr exon;
										$mpos+=abs($cdsend-$exonstart[$lastCDS_enum]);
										for (my $z=$lastCDS_enum+1;$z<$k;$z++){
											$diff+=abs($exonstart[$z]-$exonend[$z])+1;
										}
										#add the current distance
										$diff+=abs($end-$exonstart[$k]);
									}
									$mpos+=$trx_lastCDS+ $diff+1;
									$utr3{$name2}.="$name:E".sprintf("%02d",$k+1).":m.$mpos:c.+$diff;";#tested 20111111 MOOCOW
									#hv block ends
#									$utr3{$name2}++;		#positive strand for UTR3 #original
								} else {
									$exonic{$name2}.="$name:o$cdsoffset;" ;
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
								$splicing{$name2}.="$name," if ($splicing{$name2}!~/$name,/);
							}
							if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1 or $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
								$splicing{$name2}.="$name," if ($splicing{$name2}!~/$name,/);
							}
							if ($start <= $exonstart[$k] and $end>=$exonstart[$k] or $start <= $exonend[$k] and $end >= $exonend[$k]) {
								$splicing{$name2}.="$name," if ($splicing{$name2}!~/$name,/);
							}
							
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if this splicing leads to amino acid change (rather than UTR change)
							if ($splicing{$name2} and $start==$end and $start<=$cdsend) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
#									$splicing_anno{$name2} .= "$name:c.$lenexon+" . ($exonstart[$k]-$start) . revcom($ref) . '>' . revcom ($obs) . ',';
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
									my $exonnbr=$#exonstart-$k+1;
									my ($diff,$mpos);
									if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
										#query  ----
										#gene     <--*---*->
#										$utr3{$name2}++;		#negative strand for UTR5
		##hv simone
										$diff=0;$mpos=0;
										if ($start>$exonstart[$lastCDS_enum] && $start<$exonend[$lastCDS_enum]){
											#it is in the last coding exon so just calculate the distance from the cds pos
											$diff=(abs($cdsstart - $start)+1);$mpos+=abs($cdsstart-$exonend[$lastCDS_enum])+1;
										}else{
											#need to find the last coding exon UTR_length
											$diff=abs($exonstart[$lastCDS_enum]-$cdsstart)+1;#calculates the UTR of the last cds/utr exon;
											$mpos+=abs($cdsstart-$exonend[$lastCDS_enum])+1;
											for (my $z=$lastCDS_enum+1;$z<$k;$z++){
												$diff+=abs($exonstart[$z]-$exonend[$z])+1;
											}
											#add the current distance
											$diff+=abs($start-$exonend[$k]+1)+1;
										}
										$mpos+=$trx_lastCDS+ $diff;
		##
										$utr3{$name2}.="$name:E".sprintf("%02d",$exonnbr).":m.$mpos:c.+$diff;";
#										$utr3{$name2}.="$name:E".sprintf("%02d",$exonnbr).":+".abs($end-$cdsstart).";";
									} elsif ($start > $cdsend) {
										#query             ----
										#gene     <--*---*->
#										$utr5{$name2}++;		#negative strand for UTR3
										$diff=abs($start-$cdsend);
										$mpos=$cdsoffset-$diff;
										$utr5{$name2}.="$name:E".sprintf("%02d",$exonnbr).":m.$mpos:c.-$diff;";
#										$utr5{$name2}.="$name:E".sprintf("%02d",$exonnbr).":-".abs($start-$cdsend).";";
									} else {
										$exonic{$name2}.="$name:o$cdsoffset;" ;
										not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '-', $i, @exonstart-$k, $nextline];
									}
									$foundgenic++;
									last;
								} elsif ($k < @exonstart-1 and $end < $exonstart[$k+1]) { #intronic on the negative
									#HUEFLAG	
									#calculate if it is closer to the exon$k or exon$k+1
									#might need to include strand since the upstream/downstream calculator is hardcoded in the next subroutine +=downstream -=upstream but is opposite on opposite strand!!
									print "at LN".__LINE__." looking at intron..$k,$exonend[$k] and $exonstart[$k+1].\n" if ($debug);								
									my ($flg,$distToExon)=_getNearest($start,$exonend[$k],$exonstart[$k+1],'-');#flg denotes whether it's the last element or current element
									my $exonnbr=$#exonstart-$k;
									$intronic{$name2}.="$name:E".($exonnbr+$flg).":$distToExon;";#End hueflag 
									#original code									$intronic{$name2}++;
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
								my ($diff,$mpos);
								#here the trick begins to differentiate UTR versus coding exonic
								my $exonnbr=$#exonstart-$k+1;
								if ($end < $cdsstart) {			#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
									#query  ----
									#gene     <--*---*->
									$diff=0;$mpos=0;
									if ($start>=$exonstart[$lastCDS_enum] && $start<=$exonend[$lastCDS_enum]){#die "haven't coded $start\n";
										#it is in the last coding exon so just calculate the distance from the cds pos
										$diff=(abs($cdsstart - $start));$mpos+=abs($cdsstart-$exonend[$lastCDS_enum]);
									}else{
										#need to find the last coding exon UTR_length
										$diff=abs($exonstart[$lastCDS_enum]-$cdsstart)+1;#calculates the UTR of the last cds/utr exon;
										$mpos+=abs($cdsstart-$exonend[$lastCDS_enum]);
										for (my $z=$lastCDS_enum-1;$z>$k;$z--){
											$diff+=abs($exonstart[$z]-$exonend[$z]);
										}
										#add the current distance
										$diff+=abs($start-$exonend[$k]);
									}
									$mpos+=$trx_lastCDS+ $diff;
									$utr3{$name2}.="$name:E".sprintf("%02d",$exonnbr).":m.$mpos:c.+$diff;";#tested 20111111 MOOCOW  #negative strand for UTR3 tested2011/10/19
#									$utr3{$name2}.="$name:E".sprintf("%02d",$exonnbr).":+".abs($end-$cdsstart).";";
#									$utr3{$name2}++;		#negative strand for UTR5
								} elsif ($start > $cdsend) {
									#query             ----
									#gene     <--*---*->
#									$utr5{$name2}++;		#negative strand for UTR3
									$diff=abs($cdsend-$start);
									$mpos=$cdsoffset-$diff;
									$utr5{$name2}.="$name:E".sprintf("%02d",$exonnbr).":m.$mpos:c.-$diff;";
#									$utr5{$name2}.="$name:E".sprintf("%02d",$exonnbr).":-".abs($start-$cdsend).";";
								} else {
									$exonic{$name2}.="$name:o$cdsoffset;"  ;
									not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '-', $i, @exonstart-$k, $nextline];
								}
								$foundgenic=$nextgene;
								last;
							}
						}
					}
				}
			}
		}
		$foundgenic or $intergenic{''}++;		#changed $name2 to '' on 20110924
		$i =~ m/000000$/ and print STDERR "NOTICE: Finished analyzing $i query variants\n";

	
		my (@txname, %genename);
		my (%newsplicing);
		my $splice_ori = '';
		
		#process splicing annotation (change %splicing hash to %newsplicing, where gene name is replaced by gene name plus splicing annotation)
		if (%splicing) {
			if ($end-$start+1<=$splicing_threshold) {		#make sure that long indel are not considered here
				for my $tempname (keys %splicing) {
					if ($splicing_anno{$tempname}) {
						$splicing_anno{$tempname} =~ s/,$//;	#remove the trailing comma
						$tempname .= "($splicing_anno{$tempname})";
						 if (exists $gene_ori{$tempname}){$splice_ori="$gene_ori{$tempname};";}else{ $splice_ori='.';}
					}
					$newsplicing{$tempname}++;
				}
			} else {
				%newsplicing = %splicing;
			}
			
		}
		
		if ($splice_ori && $splice_ori=~/;$/){chop $splice_ori ;}else{if (keys %splicing>1){$splice_ori='.';}else{my $tmp=join(',',keys%splicing);$splice_ori=$gene_ori{$tmp};}}
		
		if ($separate) {		#separately print out each effect on one line
			if (%exonic or %splicing or %intronic or %utr5 or %utr3 or %ncrna or %upstream or %downstream) {
				my ($info,$ori);#hv made changes to the outupt
#				%exonic and ($info,$ori)=_printExtras(\%exonic) and print OUT "exonic\t", join(",", sort keys %exonic), "\t$info\t$ori\t", $nextline, "\n";
#				%newsplicing and $end-$start+1<=$splicing_threshold and print OUT "splicing\t", join (',', sort keys %newsplicing), "\t\t$splice_ori\t", $nextline, "\n";
#				%intronic and ($info,$ori)=_printExtras(\%intronic) and print OUT "intronic\t". join (",",sort keys %intronic), "\t$info\t$ori\t", $nextline, "\n";
#				%utr5 and ($info,$ori)=_printExtras(\%utr5) and print OUT "UTR5\t". join (",",sort keys %utr5)."\t$info\t$ori\t", $nextline, "\n";
#				%utr3 and ($info,$ori)=_printExtras(\%utr3) and print OUT "UTR3\t". join (",",sort keys %utr3)."\t$info\t$ori\t", $nextline, "\n";
#				%ncrna and ($info,$ori)=_printExtras(\%ncrna) and print OUT "ncRNA\t". join (",",sort keys %ncrna),"\t$info\t$ori\t$nextline\n";
#				%upstream and ($info,$ori)=_printExtras(\%upstream) and  print OUT "upstream\t". join (",",sort keys %upstream)."\t$info\t+\t", $nextline, "\n";
#				%downstream  and ($info,$ori)=_printExtras(\%downstream) and print OUT "downstream\t". join (",",sort keys %downstream)."\t$info\t+\t", $nextline, "\n";
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
							push @noncoding, $key  if (not exists $intronic{$key});
						} else {
							push @coding, $key;
						}
					}
					print EXONS "$nextline\n" if ($keepline);
#					@coding and print OUT "exonic\t", join(",", sort @coding), "\t", $nextline, "\n";
					@coding and ($info,$ori)=_printExtras(\%exonic,\@coding) and print OUT "exonic\t", join(",", sort @coding), "\t$info\t$ori\t", $nextline, "\n";
#					@noncoding and print OUT "ncRNA_exonic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
					@noncoding and ($info,$ori)=_printExtras(\%ncrna,\@noncoding) and print OUT "ncRNA_exonic\t". join(",", sort @noncoding),"\t$info\t$ori\t$nextline\n";
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
#					@coding and print OUT "splicing\t", join(",", sort @coding), "\t", $nextline, "\n";
#					@noncoding and print OUT "ncRNA_splicing\t", join(",", sort @noncoding), "\t", $nextline, "\n";
					@coding and ($info,$ori)=_printExtras(\%splicing,\@coding) and print OUT "splicing\t", join(",", sort @coding), "\t$info\t$ori\t", $nextline, "\n";
					@noncoding and ($info,$ori)=_printExtras(\%splicing,\@noncoding) and print OUT "ncRNA_splicing\t", join(",", sort @noncoding), "\t$info\t$ori\t", $nextline, "\n";
				}
				if (%intronic) {
					my (@coding, @noncoding);
					for my $key (keys %intronic) {
						if ($ncrna{$key}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
#					@coding and print OUT "intronic\t", join(",", sort @coding), "\t", $nextline, "\n";
#					@noncoding and print OUT "ncRNA_intronic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
					@coding and ($info,$ori)=_printExtras(\%intronic,\@coding) and print OUT "intronic\t", join(",", sort @coding), "\t$info\t$ori\t", $nextline, "\n";
					@noncoding and ($info,$ori)=_printExtras(\%intronic,\@noncoding) and print OUT "ncRNA_intronic\t", join(",", sort @noncoding), "\t$info\t$ori\t", $nextline, "\n";
				}
				
				#the following paragraph is commented out on 2011oct02
				#%exonic and print OUT "exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#%splicing and $end-$start+1<=$splicing_threshold and print OUT "splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#%intronic and print OUT "intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
#				%utr5 and print OUT "UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
#				%utr3 and print OUT "UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
				%utr5 and ($info,$ori)=_printExtras(\%utr5) and print OUT "UTR5\t". join (",",sort keys %utr5)."\t$info\t$ori\t", $nextline, "\n";
				%utr3 and ($info,$ori)=_printExtras(\%utr3) and print OUT "UTR3\t". join (",",sort keys %utr3)."\t$info\t$ori\t", $nextline, "\n";
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
#				%upstream and print OUT "upstream\t", join(",", sort keys %upstream), "\t", $nextline, "\n";
#				%downstream and print OUT "downstream\t", join(",", sort keys %downstream), "\t", $nextline, "\n";
				%upstream and ($info,$ori)=_printExtras(\%upstream) and  print OUT "upstream\t". join (",",sort keys %upstream)."\t$info\t+\t", $nextline, "\n";
				%downstream  and ($info,$ori)=_printExtras(\%downstream) and print OUT "downstream\t". join (",",sort keys %downstream)."\t$info\t+\t", $nextline, "\n";
			} elsif (%intergenic) {
				$genel ||= "NONE";
				$gener ||= "NONE";
				$distl ||= "NONE";
				$distr ||= "NONE";
				print OUT "intergenic\t", "$genel,$gener\t$genel(dist=$distl),$gener(dist=$distr)\t+", "\t", $nextline, "\n";
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
#				if (@coding and %splicing and $end-$start+1<=$splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
#					print OUT "exonic;splicing\t", join(",", sort @coding), ";", join (",", sort keys %newsplicing), "\t", $nextline, "\n";
#				} elsif (@coding) {
#					print OUT "exonic\t", join(",", sort @coding), "\t", $nextline, "\n";
#				} elsif (@noncoding) {
#					print OUT "ncRNA_exonic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
#				}
				print EXONS "$nextline\n" if ($keepline);
				if (@coding and %splicing and $end-$start+1<=$splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
					my ($info,$ori)=_printExtras(\%exonic,\%splicing_anno,\@coding);
					print OUT "exonic;splicing\t", join(",", sort keys %exonic), ";", join (",", sort keys %splicing_anno), "\t$info\t$ori\t", $nextline, "\n";
				} elsif (@coding) {
					my ($info,$ori)=_printExtras(\%exonic,\@coding);
					print OUT "exonic\t", join(",", sort @coding), "\t$info\t$ori\t", $nextline, "\n";
				} elsif (@noncoding){
					my ($info,$ori)=_printExtras(\%ncrna,\@noncoding);
					print OUT "ncRNA_exonic\t", join(",", sort @noncoding), "\t$info\t$ori\t", $nextline, "\n";
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
					my ($info,$ori)=_printExtras(\%splicing_anno,\@coding);
					print OUT "splicing\t", join (',', sort @coding), "\t$info\t$ori\t", $nextline, "\n";
				} elsif (@noncoding) {
					my ($info,$ori)=_printExtras(\%splicing_anno,\@noncoding);
					print OUT "splicing\t", join (',', sort @noncoding), "\t$info\t$ori\t", $nextline, "\n";
#					print OUT "ncRNA_splicing\t", join (',', sort @noncoding), "\t", $nextline, "\n";
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
				if (%utr5 or %utr3) {
					if (%utr5 and %utr3) {
						my ($info,$ori)=_printExtras(\%utr5,\%utr3,\@coding);
						print OUT "ncRNA_UTR5;ncRNA_UTR3\t", join(",",@coding), "\t$info\t$ori\t", $nextline, "\n";		#use ";" to separate UTR5 and UTR3 genes
					} elsif (%utr5) {
						my ($info,$ori)=_printExtras(\%utr5,\@coding);
						print OUT "ncRNA_UTR5\t", join(",",@coding), "\t$info\t$ori\t", $nextline, "\n";
					} else {
						my ($info,$ori)=_printExtras(\%utr3,\@coding);
						print OUT "ncRNA_UTR3\t", join(",",@coding), "\t$info\t$ori\t", $nextline, "\n";
					}
				} elsif (@noncoding) {
					my ($info,$ori)=_printExtras(\%intronic,\@noncoding);
					print OUT "ncRNA_intronic\t", join(",", sort @noncoding), "\t$info\t$ori\t", $nextline, "\n";
				} else {
					die "FATAL ERROR: please report bug to ANNOVAR author with your input line <$nextline>\n";
				}
			} elsif (%utr5 or %utr3) {
				
				if (%utr5 and %utr3) {
					my ($info,$ori)=_printExtras(\%utr5,\%utr3);
					print OUT "UTR5;UTR3\t", join(",", sort keys %utr5), ";", join(",", sort keys %utr3), "\t$info\t$ori\t", $nextline, "\n";		#use ";" to separate UTR5 and UTR3 genes
				} elsif (%utr5) {
					my ($info,$ori)=_printExtras(\%utr5);
					print OUT "UTR5\t".join (",",sort keys %utr5)."\t$info\t$ori\t", $nextline, "\n";
				} else {
					my ($info,$ori)=_printExtras(\%utr3);
					print OUT "UTR3\t".join (",",sort keys %utr3)."\t$info\t$ori\t", $nextline, "\n";
				}
			} elsif (%intronic) {
				my ($info,$ori)=_printExtras(\%intronic);
				print OUT "intronic\t" .join (",",sort keys %intronic)."\t$info\t$ori\t", $nextline, "\n";
			} elsif (%upstream or %downstream) {
				my ($info,$ori);
				if (%upstream and %downstream) {
					($info,$ori)=_printExtras(\%upstream,\%downstream);
					print OUT "upstream;downstream\t", join(",", sort keys %upstream), ";", join(",", sort keys %downstream), "\t$info\t+\t", $nextline, "\n";
				} elsif (%upstream) {
					($info,$ori)=_printExtras(\%upstream);
					print OUT "upstream\t", join(",", sort keys %upstream), "\t$info\t+\t", $nextline, "\n";
				} else {
					($info,$ori)=_printExtras(\%downstream);
					print OUT "downstream\t", join(",", sort keys %downstream), "\t$info\t+\t", $nextline, "\n";
				}
			} elsif (%intergenic) {
				$genel ||= "NONE";
				$gener ||= "NONE";
				$distl ||= "NONE";
				$distr ||= "NONE";
				print OUT "intergenic\t", "$genel,$gener", "\t$genel(dist=$distl),$gener(dist=$distr)\t+\t", $nextline, "\n";
			} else {
				die "FATAL ERROR: please report bug to ANNOVAR author with your input file\n";
			}
		}
	}
	%refseqvar and !$novar and annotateExonicVariants (\%refseqvar, $geneidmap, $cdslen, $mrnalen);#added novar 20120329
	return ($linecount, $invalidcount);
}
sub _sortFilterFile{
	my $xfn=shift;my $sfx="_filtered";
	$xfn=~s/$sfx$//g;
	eval{
		system ("cat $xfn\_filtered $xfn\_dropped > $xfn.tmp\n");
		system ("sort -n $xfn.tmp | cut -f2 > $xfn\n");
		if (!-z $xfn){
			system ("rm $xfn\_filtered $xfn\_dropped $xfn.tmp\n");
		}
	};
	if ($@){
		print "[WARN]  Could not execute _sortFile on $xfn\n";
	}
	return $xfn;
}
sub _getNearest($$$$){
	#determine which feature is nearest and determines upstream(-) or downstream (+)
	my $xquery=shift;
	my $lastElement=shift;
	my $curElement=shift;#misnomer as current element for -ve strand
	my $strand=shift;
	my $posFor=abs($curElement-$xquery);
	my $posLast=abs($xquery-$lastElement);
	print 'LN'.__LINE__."dist to last element $posLast and position to forward $posFor\n" if ($debug);
	my ($distToExon,$flg);
	if ($posFor<$posLast){#if the distance forward is greater that position back
		if ($strand eq '+'){
			$flg=0 ;#changed the flag to flg=0 for G6PC3 on 12/16/11 position chr17:42153035
			$distToExon="-$posFor";#minus is correct for positive strand
		}elsif ($strand eq '-'){
			$flg=0 ;
			$distToExon="+$posFor";}
	}else{
		if ($strand eq '+'){$flg=1; $distToExon="+$posLast";#plus sign is correct for positive strand,tested with G6PC3 12/16
		}elsif ($strand eq '-'){$flg=1 ;$distToExon="-$posLast";}
		
	}
	return ($flg,$distToExon);
}
sub _printExtras{ #additional coding information
	my $xinfo='';my $suff='';
	my @hash_arr;
	my $inp_arr;
	foreach my $i (@_){
		if (ref($i)=~/HASH/i){
			push (@hash_arr,$i);
		}elsif (ref($i)=~/ARRAY/i){
			$inp_arr=$i;
		}
	}
	if ($inp_arr){
		for (my $z=0;$z<=$#{$inp_arr};$z++){
			foreach my $href ( @hash_arr){
				while ( my ($key, $value) = each( %{$href}) ){ 
					next if ($key!~/$$inp_arr[$z]/);
					if (exists ($$href{$$inp_arr[$z]})){
						chop $value if ($value=~/[,\|:;]$/);
						$xinfo.="$value;" ;
						$suff.="$gene_ori{$key};" if (exists $gene_ori{$key} and $suff ne "$gene_ori{$key};"); 
					}
				}
			}
			
		}
	}else{
		foreach my $href (@hash_arr){
			while ( my ($key, $value) = each(%{$href}) ){ 
				chop $value if ($value=~/[,\|:;]$/);
				$xinfo.="$value;" ;
				$suff.="$gene_ori{$key};" if (exists $gene_ori{$key} and $suff ne "$gene_ori{$key};"); }
			
		}
	}
	chop $xinfo if ($xinfo=~/;$/);chop $suff if ($suff=~/;$/);
	return ($xinfo,$suff);
}
sub collapseANNOVAR{
	my $fn=shift;
	my $line=shift;
	$line||=0;
	my $currline=($hasHeader)?1:0;
	my $genes_href;
	my $exfile=$fn;$exfile=~s/variant_function/exonic_variant_function/g;
	if ($line){
		open (COLLAPSED,">>$fn.collapsed") or die "Cannot open $fn.collapsed for appending\n";
	}else{
		open (COLLAPSED,">$fn.collapsed") or die "Cannot open $fn.collapsed for writing\n";
	}
	print STDERR "printing to $fn.collapsed($line)\n";
	open (ANNOVAR,"<$fn") or die "Cannot open $fn\n";
	print STDERR "Reading $fn\n";
	my $block=0;my $href=();
	my $comment=$hasHeader;
	if ($comment ne ''){
		chomp $comment;
	}else{
		$comment='Comment';
	}
	print COLLAPSED join ("\t","#ANNOVAR annot","Gene","Rel pos to feature(s)","Gene Ori","Chr","Query Start","Query End","Allele1","Allele2","$comment")."\n" if (!$line);
	my $output='';
	LINE: while (my $line=<ANNOVAR>){
		
		next LINE if ($line>0 && $line<$currline);#if a line was specified, make sure that you read the file until it reaches that position
		if ($line=~/^#/ ){ $currline++;next;  }
		chomp $line;
		$currline++;
		if ($currline%100000==0){#write every 100000 lines
			print COLLAPSED "$output";$output='';
		}
		my @line_arr=split("\t",$line);
		push (@line_arr," ") until ($#line_arr>=9);push(@line_arr,"-1");
		if ($line=~/^exonic/ ){
			# print "working on ($currline)$line\n";<STDIN>;
			if ($currline>$block){
				$href=();#clear the hash;
				($href,$block,$genes_href)=readBlock("$exfile",$currline,$block,$genes_href);#print Dumper (\%{$href});<STDIN>;
			}
			if (!exists $$href{$currline}){
				# print Dumper ($href). "$currline";<STDIN>;
				$output.="*".join ("\t",@line_arr[0..$#line_arr-1])."\n";
			}else{
				$output.= "$line_arr[0]\t$line_arr[1]\t$$href{$currline}\t".join("\t",@line_arr[3..($#line_arr-1)])."\n";
			}
#		}elsif ($line=~/ERR/){
#			$line=~s/ERR\tERR\t/ERR\tERR\t\t/g;
		}else{
			$output.= join ("\t",@line_arr[0..$#line_arr-1])."\n";
		}
	}
	print COLLAPSED "$output" if ($output);
	close COLLAPSED;
	if ($keepline && keys (%{$genes_href})){
		open (FILE,">$fn.stats") or die "Cannot open for writing $fn.stats\n";
		my @arr=('frameshift deletion','frameshift insertion','frameshift substitution','nonframeshift deletion','nonframeshift insertion','nonframeshift substitution',
				'nonsynonymous SNV','stopgain SNV','stoploss SNV','synonymous SNV','unknown');
		print FILE join ("\t",'Gene:',@arr)."\n";
		foreach my $gene (keys %{$genes_href}){
			print FILE "$gene\t";
			foreach my $vartype (@arr){
				$$genes_href{$gene}{$vartype}||=0;
				print FILE "$$genes_href{$gene}{$vartype}\t";
				delete ($$genes_href{$gene}{$vartype});
			}
			print FILE "\n";
			if (keys %{$gene_href{$gene}}!=0){
				warn "You have the following var types to resolve for $fn and $gene\n".Dumper (\%{$gene_href{$gene}})."\n";
			}
		}
		close FILE;
	}
	return;
}
sub _collapseANNOVAR{
	my $fn=shift;
	my $currline=($hasHeader)?1:0;
	my $genes_href;
	my $exfile=$fn;$exfile=~s/variant_function/exonic_variant_function/g;
	open (COLLAPSED,">$fn.collapsed") or die "Cannot open $fn.collapsed for writing\n";
	print STDERR "printing to $fn.collapsed($line)\n";
	open (ANNOVAR,"<$fn") or die "Cannot open $fn\n";
	open (EXONIC,"<$exfile") or die "Cannot open $fn\n";
	print STDERR "Reading $fn\n";
	my $block=0;my $href=();
	my $comment=$hasHeader;
	if ($comment ne ''){
		chomp $comment;
	}else{
		$comment='Comment';
	}
	my %gene_ref;my $feat_idx;
	# if ($correctAllele){
	# 	print COLLAPSED join ("\t","#ANNOVAR annot","Gene","RefCorrect?","Rel pos to feature(s)","Gene Ori","Chr","Query Start","Query End","Allele1","Allele2","$comment")."\n";
	# 	$feat_idx=3;
	# }elsif (!$line && !$noheader){
		print COLLAPSED join ("\t","#ANNOVAR annot","Gene","Rel pos to feature(s)","Gene Ori","Chr","Query Start","Query End","Allele1","Allele2")."\n";
	# 	$feat_idx=2;
	# }
	my $output='';my $linenumb=0;my $exnumb=0;
my $debug=0;
	LINE: while (my $line=<ANNOVAR>){
		$linenumb++;
		chomp $line;
		if ($currline%100000==0){#write every 100000 lines
			print COLLAPSED "$output";$output='';
		}
		my @line_arr=split("\t",$line);
		print "reading $linenumb ($line) $fn\n" if ($debug);
		push (@line_arr," ") until ($#line_arr>=9);push(@line_arr,"-1");
		if ($line=~/^exonic/ ){
			my $exonic=<EXONIC>;chomp $exonic;$exnumb++;
			print "THis is exonic($line)\n$exonic\nContinue??($extype)" if ($debug);
			print STDERR "not matching $exonic at $fn line $linenumb and $exfile line numb $exnumb(\n\n$output\n\n)\n" and next LINE if (!$exonic);
			my @ex_arr=split("\t",$exonic);
			# $exonic=($correctAllele)?"$ex_arr[($feat_idx-1)]\t":'';
			$exonic=($extype)?"$ex_arr[1]:$ex_arr[2]":"$ex_arr[2]";
			print "\n\n$exonic??....waiting for enter" and <STDIN> if ($debug);
			# print join ("\n\t",@line_arr);<STDIN;
			print "\n\nPrinting:\n$line_arr[0]\t($line_arr[1])\t$exonic\t".join("\t",@line_arr[3..($#line_arr-2)])."\nwaiting for enter.." and <STDIN> if ($debug);
			$output.= "$line_arr[0]\t$line_arr[1]\t$exonic\t".join("\t",@line_arr[3..($#line_arr-2)])."\n";
			if ($keepline && $ex_arr[2]=~/^([\w\-\.]+):/){
				$gene_ref{$1}{$ex_arr[1]}++;
			}
#		}elsif ($line=~/ERR/){
#			$line=~s/ERR\tERR\t/ERR\tERR\t\t/g;
		}else{
			print "not exonic ($line)\nContinue??" if ($debug);
			print "\n" . join ("\t", @line_arr[0..$#line_arr-2])."\n" and <STDIN> if ($debug);
			$output.= join ("\t",@line_arr[0..$#line_arr-2])."\n";
		}
	}
	print COLLAPSED "$output" if ($output);
	close COLLAPSED;close ANNOVAR;close EXONIC;
	if ($keepline && keys (%gene_ref)){
	 	open (FILE,">$fn.stats") or die "Cannot open for writing $fn.stats\n";
	 	my @arr=('frameshift deletion','frameshift insertion','frameshift substitution','nonframeshift deletion','nonframeshift insertion','nonframeshift substitution',
	 			'nonsynonymous SNV','stopgain SNV','stoploss SNV','synonymous SNV','unknown');
	 	print FILE join ("\t",'Gene:',@arr[6..9],"nonframeshift indel/substitution", "frameshift indel/substitution",$arr[$#arr])."\n";
	 	foreach my $gene (  keys (%gene_ref)){
	 		print FILE "$gene\t";
	 		my $indel_nfs=$indel_fs=0;
	 		foreach my $vartype (@arr){
	 			$gene_ref{$gene}{$vartype}||=0;
	 			 if ($vartype=~/unknown/){
	 			 	print FILE "$indel_nfs\t$indel_fs\t$gene_ref{$gene}{$vartype}\n";
	 			 }elsif ($vartype!~/frame/i){
	 			 	print FILE "$gene_ref{$gene}{$vartype}\t"
	 			 }elsif ($vartype=~/nonframeshift/){
	 			 	$indel_nfs+=$gene_ref{$gene}{$vartype}
	 			 }else{
	 			 	$indel_fs+=$gene_ref{$gene}{$vartype};
	 			 }
	 			delete ($gene_ref{$gene}{$vartype});
	 		}
	 		
	 		if (keys %{$gene_href{$gene}}!=0){
	 			warn "You have the following var types to resolve for $fn and $gene\n".Dumper (\%{$gene_href{$gene}})."\n";
	 		}
	 	}
 	}
 	close FILE;
	 
	return;
}
sub readBlock{
	my $fn=shift;
	my $startln=shift;
	my $endln=shift;
	my $href=shift;
	my %tmphash;
	my %other_hash;
	if (!$href){
		$href=\%other_hash;
	}
	my $flg=0;
	open (FH,"<$fn") or die "cannot open $fn at LINE".__LINE__."\n";
	if ($endln<0){
		#this means I have to count the file
		$endln=~/\-(\d+)/;$endln=$1+$startln;
		$flg=1;
	}elsif ($startln>$endln){
		$endln=$startln+100;
	}
		#read in between $startln and $endln inclusive
	my $ct=0;
	while (<FH>){
		my $exln;
		my @exinfo=split("\t",$_);
		$ct++;
		if ($flg){
			if($ct >=$startln && $ct <=$endln){
				chomp;
				$tmphash{$ct}=$_ unless ($_ eq '-');
			}
		}else{
			if ($_=~/line(\d+)/){
				$exln=$1;
				if ($exln>=$startln && $exln<=$endln){
					if (exists $tmphash{$exln}){
						$tmphash{$exln}.=",$exinfo[2]";
					}else{
						if ($extype){$tmphash{$exln}="$exinfo[1]:";}
						$tmphash{$exln}.=$exinfo[2];
					}
					if ($keepline && $exinfo[2]=~/^([\w\-\.]+):/){
						$$href{$1}{$exinfo[1]}++;
					}elsif (!$keepline){
					}else{
						#warn "This shouldn't happen: (working between lines $startln && $endln) No gene on $exln and $exinfo[2]\n\t$_\n";
						#do not count
					}
				}
			}else{
				die "what is this?line(\\d+)\t[!\t]*\t([!\t]*)\t.*[\n\r]\n\t$_\n";
			}
		}
		last if ($exln>$endln);
	}
	close FH;
	return (\%tmphash,$endln,$href);
}
sub revcompl {
	my $input=shift;my $rev;
    $rev = reverse $input;
	$rev =~ tr/ACGTacgt/TGCAtgca/; 
	return $rev;
}
sub readAnnotBlock{
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
		if ($$arr[$i]=~/^#ANNOVAR/i){
			$idx{'ANVR_TYPE'}=$i;
		}elsif ($$arr[$i]=~/#(\S+)/i){
			$idx{uc($1)}=$i;
		}elsif ($$arr[$i]=~/Gene$/i){
			$idx{'ANVR_GENE'}=$i;
		}elsif ($$arr[$i]=~/Rel pos/i){
			$idx{'ANVR_ANNOT'}=$i;
		}
	}
#	print Dumper (\%idx);
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

sub smallest_substring{
	my $testseq=shift;
	my $orig=shift;
	my $shortie = (length $testseq < length($orig))?length($testseq)+1:length($orig)+1;
	my $i=1;
	if ($testseq eq $orig){
	  $i=$shortie;
	}else{
	  until (substr($orig,0,$i) ne substr($testseq,0,$i) || $i>$shortie){
	   $i++;
	  }
	}
	if ($i>1 && $i<=$shortie){
	 	return substr($orig,0,$i-1)."\n";
	}else{
		return '';
	}
}
######################## INDEXING SUBROUTINES ##########################
# usage: build_index(*DATA_HANDLE, *INDEX_HANDLE)
sub build_index {
    my $data_file  = shift;
    my $index_file = shift;
    my $block= shift;
    $block = ($block=~/^\d+/)?$block:1000;#default
    my $offset     = 0;
	my $line=0;
    while (<$data_file>) {
        print $index_file pack("N",$offset);
        $offset = tell($data_file);
    }
}

# usage: line_with_index(*DATA_HANDLE, *INDEX_HANDLE, $LINE_NUMBER)
# returns line or undef if LINE_NUMBER was out of range
sub line_with_index {
    my $data_file   = shift;	
    my $index_file  = shift;
    my $line_number = shift;

    my $size;               # size of an index entry
    my $i_offset;           # offset into the index of the entry
    my $entry;              # index entry
    my $d_offset;           # offset into the data file

    $size = length(pack("N", 0));
    $i_offset = $size * ($line_number-1);
    seek($index_file, $i_offset, 0) or return;
    print "($entry)...($size),$info\n";
    read($index_file, $entry, $size);
    $d_offset = unpack("N", $entry);
    seek($data_file, $d_offset, 0);
    return scalar(<$data_file>);
}

#hv modified to get a range of lines
sub lines_with_index {
    my $data_file   = shift;
    my $index_file  = shift;
    my $begin_line = shift;
    my $size=shift;

#    my $size;               # size of an index entry
    my $i_offset;           # offset into the index of the entry
    my $entry;              # index entry
    my $d_offset;           # offset into the data file

    $size = length(pack("N", 0));
    $i_offset = $size * ($line_number-1);
    my $info=seek($index_file, $i_offset, 1) or return $info;
    return $info;
}
#################### END OF INDEXING SUBROUTINES ####################

1;