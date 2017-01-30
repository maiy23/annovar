#!/usr/bin/perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: 497 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-11-21 13:56:59 -0800 (Mon, 21 Nov 2011) $';


our ($verbose, $help, $man);
our ($dbfile);
our ($filetype, $bin, $outfile);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'filetype=s'=>\$filetype, 'bin=i'=>\$bin, 'outfile=s'=>\$outfile) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($dbfile) = @ARGV;

$filetype ||= 'A';
$filetype =~ m/^[ABCD|(GFF)]$/ or pod2usage ("Error in argument: the -filetype argument can be only 'A' or 'B' or 'C'");
$bin ||= 100000;
$outfile ||= "$dbfile.newdb";

print STDERR "NOTICE: Two output files will be generated for use by ANNOVAR: $outfile and $outfile.idx (use -outfile to override)\n";



#step 1: generate the new output file
print STDERR "NOTICE: Running the first step of indexing (generating $outfile) ...\n";
my $foo=`/usr/bin/file $dbfile $dbfile`;chomp $foo;
if ($foo=~/CRLF/i){
	eval{
		system ("/usr/bin/dos2unix $dbfile \n");
	};
}
my $command = "echo -n > $outfile";
system ($command);
for my $i (1 .. 22, 'X', 'Y', 'M') {	
	if ($filetype eq 'A') {
		#$command = q#grep -P '^(chr)?# . $i . q#\t\d+' # . "$dbfile | sort -n -k 2 >> $outfile";
		$command = qq#grep -P '^(chr)?$i\\t\\d+' $dbfile | sort -n -k 2 >> $outfile#;
	} elsif ($filetype eq 'B') {
		$command = qq#grep -P '^\\w+\\t(chr)?$i\\t\\d+' $dbfile | sort -n -k 3 >> $outfile#;
	} elsif ($filetype eq 'C') {
		$command = qq#grep -P '^\\w+\\t\\w+\\t(chr)?$i\\t\\w+\\t\\d+' $dbfile | sort -n -k 5 >> $outfile#;
	}elsif ($filetype eq 'GFF'){#this is a file
		$command = qq#grep -P '^\\w+\\t\\w+\\t(chr)?$i\\t\\w+\\t\\d+' $dbfile | sort -n -k 4 >> $outfile#;
	}
	$verbose and print STDERR "NOTICE: Running command: $command\n";
	system ($command);
}

#step 2: generate the index file
print STDERR "NOTICE: Running the second step of indexing (generating $outfile.idx) ...\n";
$dbfile = $outfile;			#now the dbfile is the newdb generated in step 1
my $filesize = -s $dbfile;
my %region = ();
my ($offset, $lastregion, $firstoffset, $firstline) = (0, undef, 0, undef);

open (DB, $dbfile) or die "Error: cannot read from dbfile $dbfile: $!\n";
open (IDX, ">$dbfile.idx") or die "Error: cannot write to index file $dbfile.idx: $!\n";

print IDX "#BIN\t$bin\t$filesize\n";

while (my $line=<DB>) {
	my ($chr, $start);
	my $length = length ($line);
	$line=~s/[\r\n]{1,}//g;
	if ($filetype eq 'A' ) {
		($chr, $start) = split (/\t/, $line);
	} elsif ($filetype eq 'B') {
		(undef, $chr, $start) = split (/\t/, $line);
	} elsif ($filetype eq 'C') {
		(undef, undef, $chr, undef, $start) = split (/\t/, $line);
	} elsif ($filetype eq 'D'){
		(undef, undef, $chr,  $start) = split (/\t/, $line);
	}
	defined $start or die "Error: unable to find start site from input line <$line>\n";
	$start =~ m/^\d+$/ or die "Error: the start site ($start) is not a positive integer in input line <$line>\n";
	
	my $curbin = $start - ( $start % $bin );
	my $region = "$chr\t$curbin";
	$region{ $region }{ 'min' }   = $offset if (!exists $region{ $region }{ 'min' } );
	$region{ $region }{ 'max' }   = $offset + $length;
	
	$offset =~ m/000$/ and print STDERR sprintf("NOTICE: Indexing $dbfile: %d%%\r", int(100*$offset/$filesize));
	$offset += $length;
}

for my $k ( sort {$a cmp $b} keys %region ) {
	my $chr=$k;$chr=~s/chr//;
	print IDX join("\t", $chr, $region{ $k }{ 'min' }, $region{ $k }{ 'max' }), "\n";
}

print STDERR "\nDone!\n";
close(DB);
close(IDX);


=head1 SYNOPSIS

 index_annovar.pl [arguments] <db-file>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --filetype <A|B|C>		file type (default: A)
            --bin <int>			BIN size (default: 100000)
            --outfile <file>		prefix of output file name

 Function: generate index for ANNOVAR database files
 
 Example: index_annovar.pl tempdb/hg19_cg69.txt -outfile humandb/hg19_cg69.txt
          index_annovar.pl tempdb/hg19_snp131.txt -outfile humandb/hg19_snp131.txt -filetype B
 
 Version: $LastChangedDate: 2011-11-21 13:56:59 -0800 (Mon, 21 Nov 2011) $
 
 WARNING: THIS PROGRAM IS STILL IN DEVELOPMENT PHASE AND MAY CONTAIN BUGS !

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=back

=head1 DESCRIPTION

This program will generate a new database file as well as an index file, given a 
user-specified database.

The file type A, B and C are explained below:

=over 8

A: first two tab-delimited fields are chr and start

B: first three tab-delimited fields are anything, chr and start

C: first five tab-delimited fields are anything, anything, chr, anything and start

=back

=cut