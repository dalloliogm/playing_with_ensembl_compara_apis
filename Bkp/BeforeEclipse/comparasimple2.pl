#!/usr/bin/env perl
# perl script 

# Load Modules
use utf8;
use strict;
use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::Homology;
use Bio::AlignIO;


# Variabili customizzazione
my $perc_id_cutoff 	= 70; 				# Print only genes with >$perc_id_cutoff% of homology
my $default_ensembl_id 	= "ENSG00000010704";		# Default ensembl Id used
#my @ensemblIds = qw(ENSG00000010704 ENSG00000071967 ENSG00000105697 ENSG00000106327 ENSG00000110911 ENSG00000138449 ENSG00000168509);
#my @ensemblIds = qw(ENSG00000110911 ENSG00000071967);
my @ensemblIds = qw(ENSG00000110911);
#my @species = ("Homo sapiens", "Canis familiaris", "Felis catus", "Bos taurus", "Pan troglodytes", "Macaca mulatta", "Pongo pygmaeus", "Equus caballus", "Gallus gallus", "Microcebus murinus", "Myotis lucifugus", "Ochotona princeps", "Otolemur garnettii", "Tupaia belangeri");
my @species = ("Bos taurus");
#my @species = ("Homo sapiens");	# get paralogues



#Initialize Database connection
Bio::EnsEMBL::Registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org:5306/');#||die "can't connect to database";
my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Multi", 'compara','Member')||die "can't get compara adaptor\n";
my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE', $default_ensembl_id);
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');


# Lancia le funzioni principali
my %orthologues_hash = &get_all(\@ensemblIds, \@species, $perc_id_cutoff);
&write_report_file(\%orthologues_hash);

sub write_report_file
{
	# Write a tab-delimited report file with the orthologues for every gene.
	# e.g.:
	# HFE     ENSG00000010704 ENSPTRG00000017805      ENSMMUG00000023733      ENSPPYG00000016305      ENSOGAG00000007990
	# dcytb   ENSG00000071967 ENSPTRG00000012634      ENSMMUG00000007295      ENSPPYG00000012929      ENSOGAG00000000245
	#
	# Input:
	# %homologies_hash = {geneId: {specie: @(list ofensemblIds)}}

	my $orthologues_hash = shift;
#	print Dumper($orthologues_hash);
	my $report = (); 
	for my $geneId (keys %$orthologues_hash)
	{
#		print Dumper($geneId);
		print Dumper($orthologues_hash{$geneId});

#		foreach my $specie (keys %{$orthologues_hash->{$geneId}})
#		{
#		print Dumper($specie);
#			print "$specie \n";
#			# Questo é un array , non un hash!!
#		}
	}
}

sub print_alignments
{
	print;
}


sub get_all
{
	# Call &get_homologies_by_species for every ensemblId in @ensemblIds
	#
	# Inputs:
	# @$species = array of species in binomial form (e.g.: @species = ("Mus murinus", "Pan troglodytes"), in which we have to search for homologies.
	# @$ensemblIds = list of ensembl Ids
	# $perc_id_cutoff = cutoff, select only homologies with perc_id > of this value.
	#
	# Output:
	# %homologies_hash = {geneId: {specie: @(list ofensemblIds)}}
	
	my $ensemblIds = shift;
	my $species = shift;
	my $perc_id_cutoff = shift;

	my %homologies_hash = ();
	foreach my $ensemblId (@$ensemblIds)
	{
		my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE', $ensemblId);
		$homologies_hash{$ensemblId} = &get_homologies_by_specie($species, $member, $perc_id_cutoff);
	}
	return %homologies_hash;
}

sub get_homologies_by_specie
{
	# Given a Bio::EnsEMBL::Compara::Member object, a @species array, and a $perc_id_cutoff value return a list of all the homologies 
	# of the gene in the given species with a perc_id > perc_id_cutoff.
	#
	# Inputs:
	# @$species 		= array of species in binomial form (e.g.: @species = ("Mus murinus", "Pan troglodytes")
	# $member 		= B::E::Compara::Member object initialized with an ensembl Id
	# $perc_id_cutoff 	= integer, cutoff.
	#
	# Output:
	# %gene_orthologues = {$specie: @(EnsemblId)}

	my $species = shift;
	my $member = shift;
	my $perc_id_cutoff = shift;

	my %gene_orthologues = ();	# Hash which will be returned as output. %gene_orthologues = {$specie: @(EnsemblId})}

	print "\nGene: " . $member->stable_id . " (". $member->description ."\n";
	print "Showing only homologues having perc_id > $perc_id_cutoff\n";
	foreach my $specie (@$species)
	{
		my @homologies = ();
		$gene_orthologues{$specie} = @homologies;
		my $homologies = $homology_adaptor->fetch_all_by_Member_paired_species($member, $specie);

		foreach my $homology (@{$homologies})
		{

			# Print Ortholog's ensemblId and perc_id
			foreach my $member_attribute (@{$homology->get_all_Member_Attribute})
			{

				my ($member, $attribute) = @{$member_attribute};
				my $specie = $member->taxon->binomial;

				$gene_orthologues{$specie} = $homology;
#				my $perc_id = $attribute->perc_id;
				
				if ($attribute->perc_id >= $perc_id_cutoff)
				{
					unless ($specie eq "Homo sapiens")
					{
						# Aggiunge l'omologia corrente all'array $gene_orthologues{$specie} solo perc_id é suff. alto
						push(@homologies, $homology);
						print @homologies . "\n";

						print "$specie\n";
						print "  " . $member->stable_id . "\t" . $attribute->perc_id . "\t". $attribute->perc_pos ."\t" .$attribute ->perc_cov ."\n";
						$gene_orthologues{$specie} = $member->stable_id;

#						print join " ", map { $attribute->$_ } qw(perc_id perc_pos perc_cov),"\n";

#						# print alignment
#						my $simple_align = $homology->get_SimpleAlign();
#						my $alignIO = Bio::AlignIO->newFh(
#							-interleaved => 0,
#							-fh => \*STDOUT,
#							-format => "clustalw",
#							-idlength => 20);
#
#						print $alignIO $simple_align;

					}
				}
			}
		}	
#	print "\n";
	}
#	print Dumper(%gene_orthologues);
	return  %gene_orthologues;
}


sub get_all_homologies
{
	# deprecated
	# Fetch all homologies  for every specie, and then prin only those who are in @species.

	my $homologies = $homology_adaptor->fetch_all_by_Member($member);
	foreach my $homology (@{$homologies}) 
	{
		foreach my $member_attribute (@{$homology->get_all_Member_Attribute})
		{
			my ($member, $attribute) = @{$member_attribute};
			my $specie = $member->taxon->binomial;

			if (grep (/^$specie/i, @species))
			{
				print $member->taxon->binomial . ": " . $member->stable_id ."\n";
				print "\n";
			}
		}
	}
}
