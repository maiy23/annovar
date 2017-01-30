sub printTo{
	my $file=shift;
	my $info=shift;
	my $fail=0;
	open (FILE,">$file") or $fail=1;
	print FILE $info if (!$fail);
	close FILE;
	unlink($file) if ($fail);
	return $fail;
}
sub locate{
	my $exe=shift;
	if (!`which $exe`){
		if (-e "/usr/bin/$exe"){
			return "/usr/bin/$exe";
		}elsif (-e "/bin/$exe"){
			return "/bin/$exe";
		}else{
			$exe=`locate $exe | grep $exe\$ | head -n1`;chomp $exe;
		}
	}
	if ($exe){
		return $exe;
	}else{
		die "Could not run $exe from $0\n";
	}
}
sub SAFE_PROCESS{
	my $cmd=shift;
	# print "About to run $cmd\n";return 1;
	eval{
		system ("$cmd\n");
	};
	if ($?){
		return 0;
	}else{
		return 1;
	}
}

sub getOrgType{
	my $input=shift;
	my $type=shift;
	$type=lc($type);
	# my $lookup=(
	# 	'mm10'=>'ucsc','hg19'=>'ucsc','mm9'=>'ucsc',
	# 	'hg38'=>'ucsc','GRCh37'=>'refseq','GRCh38'=>'refseq',
	# 	'NCBI36'=>'refseq','NCBI37'=>'refseq','10090'=>'taxonid','9606'=>'taxonid'
	# 	);
	my %converter=(
		'mm10'=>{'ucsc'=>'mm10', 'taxonid'=>10090, 'refseq'=>'GRCh38m', 'org'=>"mouse"},
		'hg19'=>{'ucsc'=>'hg19', 'taxonid'=>9606, 'refseq'=>'GRCh37', 'org'=>"human"},
		'mm9'=>{'ucsc'=>'mm9', 'taxonid'=>10090, 'refseq'=>'NCBI36', 'org'=>"mouse"},
		'hg38'=>{'ucsc'=>'hg38', 'taxonid'=>9606, 'refseq'=>'GRCh38', 'org'=>"human"},
		'NCBI36'=>{'ucsc'=>'hg18', 'taxonid'=>9606, 'refseq'=>'NCBI36', 'org'=>"human"},
		'GRCh37'=>{'ucsc'=>'hg19', 'taxonid'=>9606, 'refseq'=>'NCBI37', 'org'=>"human"},
		'human'=>{'ucsc'=>'hg19', 'taxonid'=>9606, 'refseq'=>'GRCh37', 'org'=>"human"},
		'mouse'=>{'ucsc'=>'mm10', 'taxonid'=>10090, 'refseq'=>'GRCh38m', 'org'=>"mouse"},
	);
	if (exists $converter{$input}{$type}){
		return $converter{$input}{$type};
	}else{
		warn ("$input and $type does not exist in hash...perhaps you should add it");
		return 0;
	}
}

sub formatFasta{
	my $seq=shift;
	my $seqWithNewLines='';
	for (my $i=0;$i<=length($seq);$i+=60){
		if ( ($i+60) > length($seq)){
			$end=(length($seq)-$i);
		}else{
			$end=60;
		}
		$seqWithNewLines.=substr($seq,$i,$end). "\n";
	}
	return $seqWithNewLines;
}

sub printPBS{
	my $exe=shift;
	my $name=shift;
	my $cwd=shift;
	open (OUT ,">$cwd/PBS.$name.bat") or (record("-".__LINE__) and die "Cannot open $cwd$name.bat for writing\n");
	print OUT  "#!/bin/sh\n";
	print OUT  "#PBS -N $name \n#PBS -l cput=290:000:00,pcput=290:00:00,walltime=290:00:00\n#PBS -j oe -l pvmem=2GB,mem=2GB,ncpus=1\n#PBS_O_WORKDIR=$cwd\numask 000\n";
	print OUT "if [ ! -e \"/scratch/local/\$USER\" ];then\n\tmkdir /scratch/local/\$USER\nelse\n\tfind /scratch/local/$USER -maxdepth 1 -mtime +7 |xargs rm -rf;\nfi\n";#make temp directory if it doesn't exist on the PBS node for cdhit
	print OUT  "$exe \n";
	close OUT ;
	system ("chmod +x $cwd/PBS.$name.bat\n");
	return "$cwd/PBS.$name.bat";
}
sub runPBS{
	my $exe=shift;
	my $name=shift;
	my $cwd=shift;
	my $cmd=printPBS($exe,$name,$cwd);
	print "qsub $cmd\n";
	my $pid=`qsub $cmd`;
	return $pid;
}
1;