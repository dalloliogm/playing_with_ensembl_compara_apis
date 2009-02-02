#!/usr/bin/env perl
# perl module  by Giovanni Dall'Olio
# Module of Gene object.

package Gene;
our $AUTOLOAD;

my %fields = (
	ensemblId => undef,
	homologies => {}
);

sub new 
{
	my $classname = shift;
	my $self = {};
	
	my $self  = {
		_permitted => \%fields,
		%fields,
  	};

	bless($self, $classname);
	return $self;
}

sub AUTOLOAD 
{
	my $self = shift;
	my $type = ref($self) ||die "$self is not an object";

	my $name = $AUTOLOAD;
	$name =~ s/.*://;   # strip fully-qualified portion
	unless (exists $self->{_permitted}->{$name}) 
	{
		die "Can't access $name field in class $type";
 	}

	if (@_) {
		return $self->{$name} = shift;
 	} else {
		return $self->{$name};
 	}
}  

1;
