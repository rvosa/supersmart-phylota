package Bio::Phylo::PhyLoTA::Service::InferenceService;
use strict;
use warnings;
use Bio::Phylo::PhyLoTA::Service;
use Bio::Phylo::Util::Exceptions 'throw';
use base 'Bio::Phylo::PhyLoTA::Service';

=item create

Instantiates and configures the wrapper object, returns the instantiated wrapper.
This is an abstract method that needs to be implemented by the child class.

=cut

sub create {
    throw 'NotImplemented' => "missing create method in child class " . ref(shift);
}

=item configure

Provides the $config object to the wrapper so that settings defined in $config object can 
be applied to the wrapper. This is an abstract method that needs to be implemented by the 
child class.

=cut

sub configure {
    throw 'NotImplemented' => "missing confgure method in child class " . ref(shift);
}

=item run

Runs the analysis on the provided supermatrix. can potentially be run multiple
times (in parallel?) to do bootstrapping. returns a tree file. This is an abstract method 
that needs to be implemented by the child class.

=cut

sub run {
    throw 'NotImplemented' => "missing run method in child class " . ref(shift);
}

=item is_bayesian

Indicates whether the inference service is bayesian, in which case we will
not run a bootstrapping analysis. Note that certain wrappers, e.g. exabayes, 
need to override this

=cut

sub is_bayesian { 0 }

=item replicate

Getter/setter, is used internally to track the current bootstrap replicate

=cut

sub replicate {
	my ( $self, $rep ) = @_;
	if ( defined $rep ) {
		$self->{'replicate'} = $rep;
	}
	return $self->{'replicate'};
}


1;