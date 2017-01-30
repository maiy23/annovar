#!/usr/local/bin/perl
use Getopt::Std;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use ABCC;
use vars qw($opt_n $opt_p $opt_o $opt_c);
 getopts("p:o:n:c");
=head
This is a wrapper script to combine the two files produced by annovar.
=cut

my $usage=qq(
	$0 -p <prefix file name>
		ie. if the two ANNOVAR files are named ANNOVAR.input.variant_function and ANNOVAR.input.exonic_variant_function
		then specify \"ANNOVAR.input\" as the prefix

	[OPTIONAL]
		-o <outpufilename>
	
	[One or none of the following can be specified]  -c option takes precedent over -n option
		-n <line number to begin from> 
			if your job was not successfully completed on the first run, then specify the line in the input file where it failed
		-c <if specified, then indicates that the report was split into chunks>
			all reports with <ANNOVAR.input.variant_function_*> will be aggregated
				ANNOVAR.input.variant_function_1
				ANNOVAR.input.variant_function_2 ....
			will be combined into one file ANNOVAR.input.variant_function.collapsed (DEFAULT, or -o option ), which should not exist in the directory
);
$opt_n ||= '';
die "You must specify an input prefix\n$usage\n" if (!defined $opt_o);
if ($opt_c){
	die "haven't tested yet\n";
	$opt_o ||= "$opt_i.variant_function.collapsed";
	my @arr=split("\n",`ls $opt_i.variant_function\_*`);
	for (my $i=0;$i<=$#arr;$i++){
		collapseANNOVAR($arr[$i]);
		eval {
			system ("cat $opt_i.variant_function\_$i.collapsed >> $opt_o\n");
		};
		if ($@){
			die "[ERROR] could not concatenate $opt_i.variant_function\_$i.collapsed to $opt_o\n";
		}
	}
}else{
	collapseANNOVAR("$opt_i.variant_function");
	if (defined $opt_o){
		eval{
			system ("mv $opt_i.variant_function.collapsed $opt_o\n");
		};
		if ($@){
			die "[ERROR] There is a problem moving $opt_i.variant_function.collapsed to $opt_o\n";
		}
	}
}

print STDERR "Done!\n";
