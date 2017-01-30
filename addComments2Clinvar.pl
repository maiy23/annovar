#!/usr/bin/perl
use strict;
use Data::Dumper;
use XML::Parser;

my $bin=`dirname $0`;chomp $bin;
$/="\n\n";
my $date=`date '+%Y%m%d'`;$date=~s/[\r\n]//;
my $annovardb="/SeqIdx/annovardb/humandb/hg19_clinvar_$date.txt";
my $xml="/is1/projects/bioinfo/static/dwnld/clinvar/ClinVarFullRelease_00-latest.xml";
my $tmpDatabase = "/SeqIdx/annovardb/humandb/hg19_clinvar_$date.tmp.txt";
my $annovar_script_dir = "/bioinfoC/hue/annovar/annovar_forFDI/";
my $validate=0;#set this to print and do not consult the VCF file
my %headers = (
	'pmid'=>'PMID',
	'disease'=>"CLNDIS",
	'clinsig'=>"CLNSIG",
	'comment'=>'COMMENT',
	'clinvar_acc'=>"CLNACC"
	);
open (FILE,"<$xml") or die "Cannot open file\n";
my %hash;my $linect=0;
BLOCK: while (<FILE>){
	my $block =$_;$linect++;
	if (/<ClinVarSet/ .. /<\/ClinVarSet>/ ) {
		# my @lines=split("\n",$_);
		my ($var,$comment);
		my %matches ;#= ('disease'=>[],'pmid'=>[],'comments'=>[]);
		my @clinsigMatches;
		push (@clinsigMatches,$1) while ($_=~/(<ClinVarAssertion(.*?)(<\/ClinVarAssertion?))/gs);
		push (@clinsigMatches,$1) while ($_=~/(<ReferenceClinVarAssertion(.*?)<\/ReferenceClinVarAssertion)/gs);
		# print $1;exit;/
		# print Dumper (\@clinsigMatches);<STDIN>;
		foreach my $xblock (@clinsigMatches){
			# $xblock=~/[\s\n]{1,}/ /g;
			push (@{$matches{'disease'}},$1)  while ($xblock=~ /<TraitSet[^<>]*Type="Disease".*<Name>.*<ElementValue Type="Preferred">([^<>]*?)<\/ElementValue>[\n\s]+<\/Name>/gs);
			push (@{$matches{'disease_alt'}},$1) while ($xblock=~ /<TraitSet[^<>]*Type="Disease".*<ElementValue Type="Alternate">([^<>]*?)<\/ElementValue?/gs);
			push (@{$matches{'pmid'}},$1) while ($xblock=~ /<ID Source="PubMed">(\d+)<\/ID/gs);
			push (@{$matches{'clinvar_acc'}},"$1.$2") while ($xblock=~ /<ClinVarAccession Acc="([^"]+)" Version="(\d+)".*Type="RCV"/gs);
			push (@{$matches{'status'}},$1) while ($xblock=~ /Review Status>(.*)<\/Review/gs);
			push (@{$matches{'clinsig'}},$1) while ($xblock=~ /<\/ReviewStatus>[^<>]*<Description>([^<>]*)<\/Description/gs);
			while ($xblock=~ /<Comment>([^<>\/]*)<\/Comment/gs){
				my $f=$1;
				if ($f!~/(Converted during submission|not provided)/i){
					if($#{$matches{'comment'}}==-1){
						push(@{$matches{'comment'}},$f);
					}
				}
			}
		}
		# print Dumper (\%matches);<STDIN>;
		my @chr_loc;
		push (@chr_loc,"$1") while ($block=~/(SequenceLocation Assembly="GRCh37"(.*)?\/>)/g);
		push(@chr_loc,$1) while ($block=~/(<Attribute Change="g.\d+[-\_]\d+del\d+" Accession="NC_0+\d+" .* Type="HGVS, genomic, top level, previous" integerValue="37">)/g);
		# print "Found#of loc: $#chr_loc\n";
		
		# print "\t=>". join ("\n\t=>",@chr_loc)."|\n";<STDIN>;
		for (my $i=0;$i<=$#chr_loc;$i++){
			my $v='';
			if ($chr_loc[$i]=~/Chr="([\dXYMT]+)".*start="(\d+)" .*stop="(\d+)" .*referenceAllele="(.*)?".*alternateAllele="([^"]*)"/){
				$v=join (":",$1,$2,$3,$4,$5);
				$chr_loc[$i]=$v;
			}elsif($chr_loc[$i]=~/<Attribute Change="g.(\d+[-\_]\d+)del(\d+)" Accession="NC_0+(\d+)" .* Type="HGVS, genomic, top level, previous" integerValue="37">/){
				my $delchr=$3;my $dellen=$2;my $delloc=$1;
				my ($delstart,$delend)=split("_",$delloc);
				$v=join(":",$delchr,$delstart,$delend,"del$dellen",'-');
				$chr_loc[$i]=$v;
			}
			if ($v){
				if (exists $hash{$v}){
					foreach my $tag (keys %matches){
						push (@{$hash{$v}{$tag}},@{$matches{$tag}});
					}
				}else{
					$hash{$v}=\%matches;
				}
			}
			
		}
		#next block for debugging
		# if ($#chr_loc==-1){
		# 	print "skipping$linect;\n";
			# print $block ."\ncouldn't find sequenceLocation ";die;<STDIN>;
		# }elsif ($#chr_loc>0){
		# 	print $#chr_loc . join ("\n",@chr_loc) ."\nDumper:\n";
		# 	print Dumper ($hash{$chr_loc[0]}) . "\n";
		# 	print '================' ."\n";
		# 	print Dumper ($hash{$chr_loc[1]}) . "\n";
		# 	<STDIN>;

		# }
		# print Dumper (\%hash);<STDIN>;
	}
	# last if ($linect>100);
}

close FILE;
$/="\n";
# open (TMP,">$tmpDatabase") or die "Cannot write to tmp file at line".__LINE__."\n";
# if ($validate==0){
# 	foreach my $loc (keys %hash){
# 		my $line=join ("\t",split(":","$loc")) ."\t";
# 		foreach my $header (keys %headers){
# 			if (exists $hash{$loc}{$header}){
# 				$line.="$headers{$header}=" . uniqMe($hash{$loc}{$header}) .";";##About to print uniq vals only!!!!!
# 			}
# 		}
# 		chop $line if ($line=~/;$/);
# 		print TMP "$line\n";
# 		# print Dumper ($hash{$loc});exit;
# 	}
# }else{
# 	open (ANNOVARDB,"<$annovardb") or die "Cannot open $annovardb at line".__LINE__."\n";
# 	while (<ANNOVARDB>){
# 		next if ($_=~/^$/);
# 		my ($chr,$start,$stop,$ref,$var,$therest)=split("\t",$_);chomp $therest;
# 		my $v=join(":",$chr,$start,$stop,$ref,$var);
# 		if (exists($hash{$v}{'comment'})){
# 			print  join ("\t",$chr,$start,$stop,$ref,$var,"$therest\tCOMMENT=".join(",",@{$hash{$v}{'comment'}})) . "\n";
# 			print Dumper ($hash{$v});<STDIN>
# 		}else{
# 			print "FOO:$_";
# 			if (exists($hash{$v})){print Dumper ($hash{$v});<STDIN>;}
			
# 		}
# 	}
# 	close ANNOVARDB;
	
# 	print "Complete!  Printing results to $tmpDatabase\n";
# 	#Index and rename
	
# }
# close TMP;
my $renameDB=$tmpDatabase;$renameDB=~s/clinvar/clinvarComments/;$renameDB=~s/.tmp//;
$renameDB.="tmp2" if ($renameDB eq "$tmpDatabase");
if (-e $renameDB){
	unlink ($renameDB);
}
eval{
	system ("perl $annovar_script_dir/convert2annovarIP.pl --format bed --includeinfo  --outfile $renameDB.tmp $tmpDatabase\n");
};
if ($@){
	die "$@\n$?\n";
}
eval{
	system ("cut -f1-5,12 $renameDB.tmp >$renameDB.x\n");
	system ("perl $annovar_script_dir/index_annovar_ABCC.pl --outfile $renameDB $renameDB.x\n");
}
if ($@){
	die "$@\n$_\n";
}
if (!-z $renameDB){
	unlink ("$renameDB.tmp","x$renameDB","$tmpDatabase");
}
print "Complete!  Printing results to $renameDB and indexing\n";


sub uniqMe{
	my $val = shift(@_);
	my %hash = map {lc($_) => $_} @{$val};
	return join(",", map { "$hash{$_}" } keys %hash);
}