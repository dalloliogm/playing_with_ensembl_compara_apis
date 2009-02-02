#!/usr/bin/env perl
# perl script 

use utf8;
use strict;

use Gene;

my $gene = Gene->new();
#my $gene = Gene::new('dmt1');

$gene->ensemblId('ENSG1111');
print $gene->ensemblId();
