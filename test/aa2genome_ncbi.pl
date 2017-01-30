#!/usr/local/bin/perl
#
# translate the protein positions to genomic coordinates
# Example: KRAS:Q61L
#
#written by Uma Mudunuri
#modified for AVIA by Hue Vuong
#########################################
use Data::Dumper;
use Getopt::Std;
use vars qw($opt_h $opt_i $opt_o);
getopts("i:o:h");
my $usage=qq(
	$0 -i <input file>
[OPTIONAL]
	-h <prints this help menu
	-o <output filename>  DEFAULT STDOUT
);
print $usage and exit if ($opt_h);
my $geneFile = "/SeqIdx/annovardb/humandb/hg19_genes.gff";
open(IN5,"<$geneFile") || die ("Cannot open $geneFile\n");

my $inputFile = $opt_i;
open(IN2, $inputFile) || die ("Cannot open input file\n");
if ($opt_o){
	open (OUTPUT,">$opt_o") or die "Cannot open output file\n";
}else{
	open (OUTPUT,">&STDOUT");
}
#read mrna and exons file
while(<IN5>) {
    chomp $_;
    @line = split('\t',$_);
    $type = $line[2];
    @id = split(";", $line[8]);
    $id[0] =~ s/ID=//; $id = $id[0];

    if($id !~ /NoRefseq/) {
	$mrna = substr($id,0,index($id,".")); 
	if($type eq "mRNA") {
	    $id[1] =~ s/GeneParent=//; $rna{$id[1]} = $mrna;

	    if($orient{$mrna}) {
		#print $mrna . "\n";
	    }
	    else {
		$orient{$mrna} = $line[6]; #capture the orientation
		$chr = $line[0];
		$chr =~ s/^chr//;
		$chr{$mrna} = $chr;
	    }
	}
	if($type eq "CDS") {
	    $id =~ /_Exon(\d+)_CDS$/;
	    $exonLen = $line[4] - $line[3] + 1; 
	    push(@{$exonLen{$mrna}},$exonLen);
	    
	    $exon{$mrna}{$1}{'start'} = $line[3];
	    $exon{$mrna}{$1}{'stop'} = $line[4];
	}
    }
}
print OUTPUT "gene\tmrna\torient\tchr\tchr_start\tchr_stop\taa_start\taa_stop\n" if ($opt_h);

while(<IN2>) {
    chomp $_;
    #$line = KRAS:Q61L;
    $_ =~ /(\w+):(\w)(\d+)(\w)/;
    $featSt = $3;
    $featStop = $3;
    $gene = $1;
    $rnas = $rna{$gene};

    #print "$gene\t$rnas\t$featSt\t$featStop\n";

    @featStart = @featStop = ();
    $exonLen = $numExons = $startFlag = $stopFlag = 0; 
    if($orient{$rnas} eq "+") {
	foreach $key (sort {$a<=>$b} keys %{$exon{$rnas}}) {
	    $start = $exon{$rnas}{$key}{'start'};
	    $stop = $exon{$rnas}{$key}{'stop'};
	    $exonLen += $stop - $start + 1;
	    $prLen = int($exonLen / 3);
	    $mod = $exonLen % 3;
	    #print "$key\t$start\t$stop\t$exonLen\t$prLen\t$mod\t$featSt\t$featStop\n";
	    if(!$startFlag) {
		if(($featSt == ($prLen+1)) and ($mod > 0)) {
		    $startFlag = 1;
		    $startExon = $key;
		    $featStPos1 = ($stop - $mod) + 1;
		    $featStPos2 = $stop;
		    #print $rnas . "\n";
		}
		elsif($featSt <= $prLen) {
		    $startFlag = 1;
		    $startExon = $key;
		    $featStPos1 = $stop - (($prLen - $featSt + 1)*3) - $mod + 1;
		    $featStPos2 = $stop;				
		}
	    }
	    if($startFlag and !$stopFlag) {
		if($featStop < $prLen) {
		    $stopFlag = 1;
		    $stopExon = $key;
		    $featStopPos1 = $start;
		    $featStopPos2 = $stop - (($prLen - $featStop)*3) - $mod;
		}
	    }
	    #print "$featStPos1\t$featStPos2\t$featStopPos1\t$featStopPos2\n";
	}
	if($startExon == $stopExon) {
	    $featStart[$numExons] = $featStPos1;
	    $featStop[$numExons] = $featStopPos2;
	}
	else {
	    $featStart[$numExons] = $featStPos1;
	    $featStop[$numExons] = $featStPos2;
	    for($i=$startExon+1; $i<$stopExon; $i++) {
		$numExons++;
		$featStart[$numExons] = $exon{$rnas}{$i}{'start'};
		$featStop[$numExons] = $exon{$rnas}{$i}{'stop'};
	    }
	    $numExons++;
	    $featStart[$numExons] = $featStopPos1;
	    $featStop[$numExons] = $featStopPos2;
	}
    }
    elsif ($orient{$rnas} eq "-") {
	foreach $key (sort {$a<=>$b} keys %{$exon{$rnas}}) {
	    $stop = $exon{$rnas}{$key}{'start'};
	    $start = $exon{$rnas}{$key}{'stop'};
	    $exonLen += $start - $stop + 1;
	    $prLen = int($exonLen / 3);
	    $mod = $exonLen % 3;
	    #print "$key\t$start\t$stop\t$exonLen\t$prLen\t$mod\t$featSt\t$featStop\n";
	    if(!$startFlag) {
		if(($featSt == ($prLen+1)) and ($mod > 0)) {
		    $startFlag = 1;
		    $startExon = $key;
		    $featStPos1 = $stop;
		    $featStPos2 = $stop + $mod - 1;
		    #print $rnas . "\n";
		}
		elsif($featSt <= $prLen) {
		    $startFlag = 1;
		    $startExon = $key;
		    $featStPos1 = $stop;
		    $featStPos2 = $stop + (($prLen - $featSt + 1)*3) + $mod - 1;				
		}
	    }
	    if($startFlag and !$stopFlag) {
		if($featStop < $prLen) {
		    $stopFlag = 1;
		    $stopExon = $key;
		    $featStopPos1 = $stop + (($prLen - $featStop)*3) + $mod;
		    $featStopPos2 = $start;
		}
	    }
	    #print "$featStPos1\t$featStPos2\t$featStopPos1\t$featStopPos2\n";
	}
	if($startExon == $stopExon) {
	    $featStart[$numExons] = $featStopPos1;
	    $featStop[$numExons] = $featStPos2;
	}
	else {
	    $featStart[$numExons] = $featStPos1;
	    $featStop[$numExons] = $featStPos2;
	    for($i=$startExon+1; $i<$stopExon; $i++) {
		$numExons++;
		$featStart[$numExons] = $exon{$rnas}{$i}{'start'};
		$featStop[$numExons] = $exon{$rnas}{$i}{'stop'};
	    }
	    $numExons++;
	    $featStart[$numExons] = $featStopPos1;
	    $featStop[$numExons] = $featStopPos2;
	}
    }
    else {
	print "Error for $gene - $rnas: no orientation\n";
    }
    for($j=0; $j<=$#featStart; $j++) {
	$start = $featStart[$j];
	$stop = $featStop[$j];
	if ($opt_h){
		print OUTPUT "$gene\t$rnas\t$orient{$rnas}\t$chr{$rnas}\t$start\t$stop\t$featSt\t$featStop\n"; 
    }else{
    	print OUTPUT "$chr{$rnas}\t$start\t$stop\t.\t.\t$_\n"; 
    }
}
}
exit 0;

