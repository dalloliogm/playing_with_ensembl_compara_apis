#!/usr/bin/env perl
# Dato un ensemblId e la specie corrispondente, scarica la sequenza tramite le API di ensembl.

# moduli standard
use utf8;
use strict;
use warnings;

# librerie specifiche
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use Bio::Seq;
use Set::Scalar;
use Getopt::Long;


my $usage = qq{
	perl get_sequences.pl 

	--idfile 
		file containing ensembl Ids to retrieve.	
}


# Leggi argomenti da shell
my $idfile = shift||'./Data/Homologies/homologies.txt';
my $FastaDir = shift||'./Data/Sequences/'; 			# Where to output fasta sequences

# dichiarazione variabili 
my $registry = 'Bio::EnsEMBL::Registry';
#my $idfile = 'ids.txt'; 				# File che contiene la lista degli id, per adesso.
#my @species = ('Human', 'Chimp', 'Rhesus', 'Orangutan', 'Bushbaby');
#my %IdHash = ();


# ----- #


# Leggi il file con la tabella di omologhi e stampalo su schermo
my ($IdHash, $species) = &readIdFile;
my %IdHash = %$IdHash;
#print $IdHash->{'ENSG00000110911'}->{'horse'};
#print "IDcow: " . $IdHash{'ENSG00000110911'}{'Bos taurus'}."\n";
#map {print "$_\t"} @{$IdHash{'ENSG00000110911'}{'Bos taurus'}};
print "\n";
&printIdHash(\%IdHash);

# Connettiti al database ensembl e scrivi le sequenze su file
my %Genes = &dbconnect;
&writeFastas(\%Genes);


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
	# Returns %Genes, an hash of Bio::EnsEMBL::Gene object for every gene and for every specie.
	# %Genes: {GeneId: {Specie: <Bio::EnsEMBL::Gene> Object}}
	
	$registry->load_registry_from_db(
	    -host => 'ensembldb.ensembl.org',
	    -user => 'anonymous'
	);
	my %gene_adaptors = ();
	my %protein_adaptors = ();
	while (defined(my $organism = $species->each))
	{
		$gene_adaptors{$organism} = $registry->get_adaptor($organism, 'Core', 'Gene');
		$protein_adaptors{$organism} = $registry->get_adaptor($organism, 'Core', 'Protein');
	}

	my %Genes = ();
#	print Dumper($IdHash);
	foreach my $geneId (keys (%IdHash))
	{
		$Genes{$geneId} = ();
		print "Getting slices for gene: $geneId \n";
#		foreach my $organism ($species)
		while (defined(my $organism = $species->each))		# Cycling on a set. see perldoc Set::Static.
		{
#			my $ensemblId = $IdHash{$geneId}{$organism};
			my $IdList = $IdHash{$geneId}{$organism};
#			print "organism: $organism\n";
#			print "organism_value: " . $IdHash{$geneId}{$organism} ."\n";
			foreach my $ensemblId (@{$IdHash{$geneId}{$organism}})
			{
#				print "ensemblId: $ensemblId\n";
				unless ($ensemblId =~ m/Undef/)
				{
					my $gene_slice = $gene_adaptors{$organism}->fetch_by_stable_id($ensemblId);
					if (exists $Genes{$geneId}{$organism})
					{
						push @{$Genes{$geneId}{$organism}}, $gene_slice;
					}
					else
					{
						$Genes{$geneId}{$organism} = [$gene_slice];
					}
#					my $seq = $gene_slice->seq();
#					$seq = substr($seq, 0, 80);
#					print "$seq\n",
				}
			}
		}
#		print "\n";
	}
	return %Genes;
}

sub writeFastas
{
	# For every gene, write a Fasta file with the sequences of every homolog.
	# input:
	# \%Genes hash created by &dbconnect
	
	my $Genes = shift;
#	print $Genes;

	foreach my $geneId (keys %$Genes)
	{
#		my $FastaDir = './Data/Sequences/';
		my $genefastafile = '>' . $FastaDir . "Genes/" . $geneId . ".fasta";
		my $genefastafile = '>' . $FastaDir . "Proteins/" . $geneId . ".fasta";

		my $geneFastaIO = Bio::SeqIO->new(-file => $genefastafile, -format => 'fasta');
		my $proteinFastaIO = Bio::SeqIO->new(-file => $genefastafile, -format => 'fasta');

		for my $specie (keys %{$Genes->{$geneId}})
		{
#			print "specie: $specie\n";
#			print "specie_value: ".$Genes{$geneId}{$specie}."\n";
			foreach my $gene_slice (@{$Genes{$geneId}{$specie}})
			{
#				my $gene_slice = $Genes->{$geneId}->{$specie};
#				print Dumper($gene_slice);
				my $SeqObj = Bio::Seq->new(-seq => $gene_slice->seq(), -display_id => $gene_slice->display_id);
				print "writing sequence of $geneId (".$gene_slice->stable_id.") from $specie in >$genefastafile\n";
				$geneFastaIO->write_seq($SeqObj);
				$proteinFastaIO->write_seq($SeqObj);
			}
		}
	}
}

sub readIdFile{
	# Read the file ids.txt and create the array with the IDs. 
	#
	# Output:
	# %IdHash = {dcytb: {'Human': 'ENSG00000071967', 'Chimp': 'ENSPTRG00000017805'...}}

	my %IdHash;# = ();
	my $species = new Set::Scalar();

	open IDFILE, "< $idfile";
	while (<IDFILE>)
	{
		unless ($_ =~ m/^[ #]/)
			{
#			print $_ . "\n";
#			my @species = @species;		# reference outer var @species with an internal var @species
			chomp;
			my @fields = split(/\t/, $_);

#			my $display_label = shift(@fields);
			my $geneId = shift(@fields);
#			print "\n$geneId\n";
			$IdHash{$geneId} = {};

			foreach my $field (@fields)
				{
					print "field: $field\n";
					my ($current_specie, $id, $protein_id) = split(/:/, $field);
					$current_specie =~ s/\"//g; 
					$species->insert($current_specie);
					if (exists $IdHash{$geneId}{$current_specie})
					{
						push (@{$IdHash{$geneId}->{$current_specie}}, $id);
#					print $id, "\t", $current_specie, "\n";
					}
					else
					{
						$IdHash{$geneId}->{$current_specie} = [$id];
#					print $current_specie, "\t";
#					print $id, "\n";
					}
				}
			}
		else 
		{
			print "commento o linea vuota: $_\n";
		}

	}
	close IDFILE;
	$species->delete('""');
#	print $species;
#	print Dumper(%IdHash);
	return (\%IdHash, $species);
}

sub printIdHash
{
	# Pretty print of %IdHash created by &readIdFile.

	my $IdHash = shift;
#	my %IdHash = %$IdHash;

	foreach my $geneId (sort keys %$IdHash)
	{
		print "geneId: $geneId \n";
		for my $specie (keys %{$IdHash->{$geneId}})
		{
			print "  $specie: ";# $IdHash->{$geneId}->{$specie}\n";
			foreach my $ensemblId (@{$IdHash->{$geneId}->{$specie}})
			{
				print " $ensemblId\t";
			}
			print "\n";
		}
		print "\n";
	}
}
