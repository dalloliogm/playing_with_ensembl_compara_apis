#!/usr/bin/env perl
# Dato un ensemblId e la specie corrispondente, scarica la sequenza tramite le API di ensembl.

# moduli standard
use utf8;
use strict;
use warnings;

# librerie specifiche
use Bio::EnsEMBL::Registry;

# dichiarazione variabili 
my $registry = 'Bio::EnsEMBL::Registry';
my $idfile = 'ids.txt'; 				# File che contiene la lista degli id, per adesso.
my @species = ('Human', 'Chimp', 'Rhesus', 'Orangutan', 'Bushbaby');
my %IdHash = ();
my $FastaDir = '../Data/Sequences/Genes'; 			# Where to output fasta sequences


# ----- #


# Leggi il file con la tabella di omologhi e stampalo su schermo
my %IdHash = &getIdHash;
&printIdHash(\%IdHash);

# Connettiti al database ensembl e scrivi le sequenze su file
my %Slices = &dbconnect;
&writeFastas(\%Slices);


sub feature2string
{
	# pretty print of a slice/gene
	my $feature = shift;

	my $stable_id  = $feature->stable_id();
	my $seq_region = $feature->slice->seq_region_name();
	my $start      = $feature->start();
	my $end        = $feature->end();
	my $strand     = $feature->strand();

	return sprintf( "%s: %s:%d-%d (%+d)", $stable_id, $seq_region, $start, $end, $strand );
}

# connetti ad ensembl
sub dbconnect
{
	# Connect to ensembl database and get an adaptor for every specie.
	# Returns %Slices, an hash of Bio::EnsEMBL::Gene object for every gene and for every specie.
	# %Slices: {GeneId: {Specie: <Bio::EnsEMBL::Gene> Object}}
	
	$registry->load_registry_from_db(
	    -host => 'ensembldb.ensembl.org',
	    -user => 'anonymous'
	);
	my %adaptors = ();
	foreach my $specie (@species)
	{
		$adaptors{$specie} = $registry->get_adaptor($specie, 'Core', 'Gene');
	}

	my %Slices = ();
	foreach my $geneId (keys (%IdHash))
	{
		$Slices{$geneId} = ();
		print "Getting slices for gene: $geneId \n";
		foreach my $organism (@species)
		{
			my $ensemblId = $IdHash{$geneId}{$organism};
			unless ($ensemblId =~ m/Undef/)
			{
				my $query_slice = $adaptors{$organism}->fetch_by_stable_id($ensemblId);
				$Slices{$geneId}{$organism} = $query_slice;
				my $seq = $query_slice->seq();
				$seq = substr($seq, 0, 80);
			}
		}
		print "\n";
	}
	return %Slices;
}

sub writeFastas
{
	# For every gene, write a Fasta file with the sequences of every homologue.
	# input:
	# \%Slices hash created by &dbconnect
	
	my $Slices = shift;
#	print $Slices;
	foreach my $geneId (keys %$Slices)
	{
		my $filename = "./Sequences/$geneId.fasta";
		print "writing $filename\n";
		open FASTAFILE, ">", $filename;

		for my $specie (keys %{$Slices->{$geneId}})
		{
			my $slice = $Slices->{$geneId}->{$specie};
#			print $slice->stable_id ."\n";
#			print $slice->source ."\n";
#			print $slice->specie ."\n";
			my $header = ">" . $slice->stable_id . ";" . $specie . "\n";
			my $seq = $slice->seq() . "\n";
			print FASTAFILE $header;
			print FASTAFILE $seq;
		}
	}
}

sub getIdHash{
	# Read the file ids.txt and create the array with the IDs. 
	# %IdHash = {dcytb: {'Human': 'ENSG00000071967', 'Chimp': 'ENSPTRG00000017805'...}}
	my %IdHash;# = ();

	open IDFILE, "< $idfile";
	while (<IDFILE>)
	{
		unless ($_ =~ m/^[ #]/)
			{
#			print $_ . "\n";
			my @species = @species;		# reference outer var @species with an internal var @species
			my @fields = split();

			my $geneId = shift(@fields);
#			print "\n$geneId\n";

			foreach my $id (@fields)
				{
				my $current_specie = shift(@species);
				$IdHash{$geneId}{$current_specie} = $id;
#				print "  $current_specie: $id \n";  
				}
		print $IdHash{'HFE'}{"Human"} eq 'ENSG00000010704'|| die "Check datafile\n";
			}
		else 
		{
			print "commento o linea vuota: $_\n";
		}

	}
	close IDFILE;
#	print %IdHash;
	return %IdHash;
}

sub printIdHash
	{
		# Pretty print of %IdHash created by &getIdHash.

		my $IdHash = shift;
#		my %IdHash = %$IdHash;

		foreach my $geneId (sort keys %$IdHash)
		{
			print "geneId: $geneId \n";
			for my $specie (keys %{$IdHash->{$geneId}})
			{
				print "  $specie: $IdHash->{$geneId}->{$specie}\n";
			}
			print "\n";
		}
	}
