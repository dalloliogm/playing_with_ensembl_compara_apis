#!/usr/bin/env perl
# perl script 

use utf8;
use strict;
use warnings; 
use Data::Dumper;
use Bio::AlignIO;
use Getopt::Long;
use Path::Abstract;

# Read arguments
my $usage = qq{
	perl convert2nexus.pl
	Converts a list of clustalw files to nexus format.

		Getting help:
			[--help]

		inputfiles:
			[--inputfiles inputfiles]
			      	input files. default is ./Tests/Data/HAMP.aln
		outputfiledir:
			[--outputfiledir outputfiledir]
				directory where outputfiles will be saved, in nexus format. 
				Files will keep the same original name but their extension will be changed into .nexus.

};

my $help = undef;
my @inputfiles = ();#'./Tests/Data/HAMP.aln';
my $outputfiledir = './Tests/Data/HAMP.nexus';
my $inputformat = "clustalw";
my $outputformat = "nexus";

GetOptions(
	"help" => \$help,
	"inputfiles=s{1,}" => \@inputfiles,
	"outputfiledir=s" => \$outputfiledir,
	"inputformat:s" => \$inputformat,
	"outputformat:s" => \$outputformat
);

if ($help)
{
	print $usage;
	exit(0);
};

&read_align;

sub read_align
{
	foreach my $inputfilename (@inputfiles)
	{
		my $outputfilename = Path::Abstract->new($inputfilename)->last;
		$outputfilename =~ s/\..*$//;
		$outputfilename = $outputfiledir . "/" . $outputfilename . ".$outputformat";
#		print $outputfilenae . "\n";	

		print "reading $inputformat file $inputfilename\n";
		my $in = Bio::AlignIO->new(-file => "<".$inputfilename, -format => $inputformat, -show_symbols => 0, -show_endblock => 0)||die "can't read $inputfilename file\n";
		my $out = Bio::AlignIO->new(-file => ">".$outputfilename, -format => $outputformat, -show_symbols => 0, -show_endblock => 0)||die "can't write $outputfilename file\n";

		print "writing $outputformat file $outputfilename\n";
		while ( my $aln = $in->next_aln() )
		{
			$out->write_aln($aln)
		}
	}
}
