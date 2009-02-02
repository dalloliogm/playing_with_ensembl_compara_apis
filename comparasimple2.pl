#/usr/bin/env perl
# Given a list of ensemblIds, get their homologies from compara and write a report file on $reportfiledir. 

# Load Modules
use utf8;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::Homology;
use Bio::AlignIO;
use Bio::SeqIO;

# Opzioni dalla riga di comando (in genere Make)
my $usage = qq{
	perl comparasimple.pl 
	Given a list of ensemblIds, get their homologues from the compara database at ensembl and write a report file on \$repordir

		species:
			species in which you want to find the orthologs. The reference specie is always H. sapiens.
			Use binomial names delimited by "" or ''' (e.g. "Homo sapiens" "Canis familiaris").

		ensemblIds:
			ensembl ids of protein or transcript of which you want to find the orthologs.

		[ensemblIdtype]: (not implemented)
			Type of ensembl ID to retrieve. Could be ENSEMBLGENE, ENSEMBLTRAN or ENSEMBLPEP. default ENSEMBLGENE.

		[reportfiledir]:
			directory where report files will be saved. Default ./Tests/Data.

		[alignmentsdir]:
			If defined (default undef) write the alignments for every hortholog pair in alignmentsdir.

		[perc_id_cutoff]:
			Default percent Id to determine whether an hortholog is a false positive or not. Default 70.
};

my @species = ();
my @ensemblIds = ();
my $reportfiledir = "./Tests/Data/";
my $alignmentsdir = undef;
my $perc_id_cutoff = 70; 				# Print only genes with >$perc_id_cutoff% of homology
my $ensemblIdtype = "ENSEMBLGENE";

GetOptions(
	"help" => \$usage,
	"reportfiledir=s" => \$reportfiledir,
	"alignmentsdir:s" => \$alignmentsdir,
	"perc_id_cutoff:i" => \$perc_id_cutoff,
	"ensemblIds=s{1,}" => \@ensemblIds,
	"ensemblIdtype:s" => \$ensemblIdtype,
	"species=s{1,}" => \@species,
	);

#print "@species -\n";
#print "@ensemblIds\n";
#print "$ensemblIdtype\n";

#{
#	@species = ("Homo sapiens", "Canis familiaris", "Felis catus", "Bos taurus", "Pan troglodytes", "Macaca mulatta", "Pongo pygmaeus", "Equus caballus", "Gallus gallus", "Microcebus murinus", "Myotis lucifugus", "Ochotona princeps", "Otolemur garnettii", "Tupaia belangeri");
#	@ensemblIds = qw(ENSG00000010704 ENSG00000071967 ENSG00000105697 ENSG00000106327 ENSG00000110911 ENSG00000138449 ENSG00000168509);
#}

#Initialize Database connection
Bio::EnsEMBL::Registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org:5306/');#||print "can't connect to database";
my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Multi", 'compara','Member')||die "can't get compara adaptor\n";
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');


# Lancia le funzioni principali
my %homologies_hash = &get_all(\@ensemblIds, \@species, $perc_id_cutoff, $reportfiledir);
&write_report(\%homologies_hash, $reportfiledir);
if (defined($alignmentsdir))
{
	&write_alignments(\%homologies_hash);
}

sub write_alignments
{
	# Write the alignments of every gene with its candidate orthologs.
	# For every gene, a file '../Data/Alignments/EnsemblOrthologs/$ensemblId.aln' is written.
	#
	# Input:
	# %homologies_hash = {geneId: {specie: @(list of homology objects)}}

	print "Writing pair-wise alignments\n";

	my $homologies_hash = shift;
#	print Dumper($homologies_hash);

	for my $geneId (keys %$homologies_hash)
	{
		
		my $outputfile = "> $alignmentsdir/$geneId.aln"; 
		print "outputfile: $outputfile\n";
		open ALIGNMENTFILE, $outputfile;

		# Use a bioperl Bio::AlignIO module to handle alignment format
		my $alignIO = Bio::AlignIO->newFh(
		    -interleaved => 0,
		    -fh => \*ALIGNMENTFILE,
		    -format => "clustalw",
		    -idlength => 20);

		foreach my $specie (keys %{$homologies_hash->{$geneId}})
		{
			print ALIGNMENTFILE "\n################### $specie ##################\n\n";
			my $homologies = $homologies_hash{$geneId}{$specie};
#			print "homologies: @homologies \n";
#			print " @homologies ";
			foreach my $homology (@$homologies)
			{
#				print "  id: ". $homology->description ."\n";
				my $simple_align = $homology->get_SimpleAlign();
				print $alignIO $simple_align;

			}
			print "\n";
		}
		close ALIGNMENTFILE;
	}
}

sub write_report
{
	# Write the fasta sequence for every EnsemblID 
	# For every gene, a file '../Data/Sequences/ByEnsemblId/$ensemblId.fasta' is written.
	#
	# Input:
	# %homologies_hash = {geneId: {specie: @(list of homology objects)}}

	my $homologies_hash = shift;
	my $reportfiledir = shift;

	print "Writing report file on $reportfiledir/homologies_report.txt";

	# Open report file handle
	my $outputfile = "> $reportfiledir/homologies_report.txt";
	open REPORTFILE, $outputfile;
#	open REPORTFILE, '>&STDOUT';

	for my $geneId (keys %$homologies_hash)
	{
		my $firstcolumn_isprinted = 0;

		foreach my $specie (keys %{$homologies_hash->{$geneId}})
		{
			my $homologies = $homologies_hash{$geneId}{$specie};
#			print "homologies: @homologies \n";
#			print " @homologies ";
			foreach my $homology (@$homologies)
			{
#				print Dumper($homology);
				foreach my $member_attribute (@{$homology->get_all_Member_Attribute})
				{
					my ($member, $attribute) = @{$member_attribute};
					my $protein_member = $member->get_longest_peptide_Member;

					if ($firstcolumn_isprinted == 0)
					{
						$protein_member = $member->get_longest_peptide_Member;
						print REPORTFILE $member->display_label ."\t";
#						print $protein_member->stable_id;
						print REPORTFILE "\"Homo sapiens\":" . $member->stable_id . ":" . $protein_member->stable_id . "\t";
#						print REPORTFILE '"Homo sapiens"' . ":" . $geneId . "\t";
						$firstcolumn_isprinted = 1;
					}

#					print $member->display_label;
					unless ($member->taxon->binomial eq 'Homo sapiens')
					{
						print REPORTFILE "\"".$member->taxon->binomial."\"". ":" .$member->stable_id.":".$protein_member->stable_id."\t";
#						print REPORTFILE "\"".$protein_member->taxon->binomial."\"". ":" .$protein_member->stable_id."\t";
#						print Dumper($member);
					}
				}
			}
		}
		print REPORTFILE "\n";
	}
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
	my $reportfiledir = shift;

	my %homologies_hash = ();
	foreach my $ensemblId (@$ensemblIds)
	{
#		my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE', $ensemblId);
		my $member = $member_adaptor->fetch_by_source_stable_id($ensemblIdtype, $ensemblId);
		print $member;
		$homologies_hash{$ensemblId} = &get_homologies_by_specie($species, $member, $perc_id_cutoff, $reportfiledir);
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
	# %gene_orthologs = {$specie: @(Homology object)}

	my $species = shift;
	my $member = shift;
	my $perc_id_cutoff = shift;
	my $reportfiledir = shift;
	open REPORTFILEWITHPERCENT, ">$reportfiledir/homologies_with_percent.txt";

	my %gene_orthologs = ();	# Hash which will be returned as output. %gene_orthologs = {$specie: @(EnsemblId})}

	print REPORTFILEWITHPERCENT "\nGene: " . $member->stable_id . " (". $member->display_label .")\n";
#	print "Showing only homologues having perc_id > $perc_id_cutoff\n";
	foreach my $specie (@$species)
	{
		my @homologies = ();
		$gene_orthologs{$specie} = ();
		my $homologies = $homology_adaptor->fetch_all_by_Member_paired_species($member, $specie);

		foreach my $homology (@{$homologies})
		{
			# Print Ortholog's ensemblId and perc_id
			foreach my $member_attribute (@{$homology->get_all_Member_Attribute})
			{
				my ($member, $attribute) = @{$member_attribute};
				my $specie = $member->taxon->binomial;
				my $ensemblId = $member->stable_id;

#				my $perc_id = $attribute->perc_id;
				
				if ($attribute->perc_id >= $perc_id_cutoff)
				{
					unless ($specie eq "Homo sapiens")
					{
						# Aggiunge l'omologia corrente all'array $gene_orthologs{$specie} solo perc_id é suff. alto
#						push(@homologies, $ensemblId);
						push(@homologies, $homology);
						$gene_orthologs{$specie} = \@homologies;

						print REPORTFILEWITHPERCENT "  " . $ensemblId . "\t" . $attribute->perc_id . "\t". $attribute->perc_pos ."\t" .$attribute ->perc_cov ."\t";
						print REPORTFILEWITHPERCENT "$specie\n";
					}
				}
			}
		}	
#		print " @homologies ";
#	print "\n";
	}
#	print Dumper(%gene_orthologs);
	return  \%gene_orthologs;
}


