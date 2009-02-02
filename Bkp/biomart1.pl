
# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use strict;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my $confFile = "http://www.biomart.org/biomart/martservice?type=registry";
#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#

my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

		
	$query->setDataset("hsapiens_gene_ensembl");
	$query->addFilter("ensembl_gene_id", ["ENSG00000010704","ENSG00000138449","ENSG00000110911","ENSG00000168509","ENSG00000071967","ENSG00000106327"]);
	$query->addAttribute("ensembl_gene_id");
	$query->addAttribute("human_paralog_ensembl_gene");
	$query->addAttribute("human_paralog_percent_coverage");
	$query->addAttribute("human_paralog_percent_identity");
	$query->addAttribute("cat_orthology_type");
	$query->addAttribute("cat_percent_identity");
	$query->addAttribute("cat_homolog_percent_identity");
	$query->addAttribute("chimp_ensembl_gene");
	$query->addAttribute("chimp_percent_identity");
	$query->addAttribute("chimp_homolog_percent_identity");
	$query->addAttribute("cow_ensembl_gene");
	$query->addAttribute("cow_percent_identity");
	$query->addAttribute("cow_homolog_percent_identity");

$query->formatter("TSV");

my $query_runner = BioMart::QueryRunner->new();
############################## GET COUNT ############################
# $query->count(1);
# $query_runner->execute($query);
# print $query_runner->getCount();
#####################################################################


############################## GET RESULTS ##########################
# to obtain unique rows only
# $query_runner->uniqueRowsOnly(1);

$query_runner->execute($query);
$query_runner->printHeader();
$query_runner->printResults();
$query_runner->printFooter();
#####################################################################
