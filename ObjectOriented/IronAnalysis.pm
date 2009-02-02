#!/usr/bin/env perl
# perl module  by Giovanni Dall'Olio

package IronAnalysis;

sub new
{
	my $caller = shift;
	my $class = ref($caller) || $caller;
	  
	my $this = {};

	$this->{'genes'} = {}	# Array of Gene objects

	bless ($this, $class);
	return $this;
}

sub getGeneByName
{
	my $this = shift;
	my $gene_name = shift;

	print ;
}

sub getGeneByEnsemblId
{
	my $this = shift;
	my $ensemblId = shift;

	print ;
}

sub addGeneObject
{
	print ;
}

return 1;

