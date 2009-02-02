#!/usr/bin/env perl
# Given an EnsemblId, get its fasta sequence from ensembl
# Usage:
# getFastaSequenceFromEnsembl ENSBTAP00000003041 Cow -p

use utf8;
use strict;
use warnings;
use Data::Dumper;

my $ensemblId = $ARGV[0];# ||$ensemblId = "ENSBTAP00000003041";
my $specie = $ARGV[1];# ||$ensemblId = "Bos taurus";

#print "ensemblId: $ensemblId \n";
#print "specie: $specie\n";

use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous'
	);#|| die "can't load registry";
#$registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org:5306/')||die 'cant load registry\n';

my $gene_adaptor = $registry->get_adaptor($specie, 'Core', 'Gene');
#$gene_adaptor->dbc();
#print Dumper($gene_adaptor);

my $gene = $gene_adaptor->fetch_by_stable_id($ensemblId);
#print Dumper($gene);
#print $gene;

#print $gene->common_name;
#print $gene->seq();
print ">" . $gene->stable_id . "\n";
print $gene->seq();
print "\n";

