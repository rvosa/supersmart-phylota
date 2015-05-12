package Bio::Phylo::PhyLoTA::Service;
use strict;
use warnings;
use Moose;
use Bio::Phylo::PhyLoTA::DAO;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::Logger;
use Bio::Phylo::Util::Exceptions 'throw';

my $config = Bio::Phylo::PhyLoTA::Config->new;
my $schema = Bio::Phylo::PhyLoTA::DAO->new;
my $logger = Bio::Phylo::Util::Logger->new;

=head1 NAME

Bio::Phylo::PhyLoTA::Service - base class for the service layer

=head1 DESCRIPTION

The functionality of the pipeline is primarily implemented by the service layer, i.e.
by packages inside this namespace. The package will be subdivided by domain, e.g.
having to do with sequences, taxa, trees, and fossils, but they will share functionality,
which is implemented here.

=head1 METHODS

=over

=item schema

Returns the L<Bio::Phylo::PhyLoTA::DAO> singleton

=cut

sub schema { $schema }

=item config

Returns the L<Bio::Phylo::PhyLoTA::Config> singleton

=cut

sub config { $config }

=item logger

Returns the L<Bio::Phylo::Util::Logger> singleton

=cut

sub logger { $logger }

=item find_seq

Given a sequence GI, returns the L<Bio::Phylo::PhyLoTA::DAO::Result::Seq> object.

=cut

sub find_seq {
	my ( $self, $gi ) = @_;
	my $result;
	eval {	        	        
	        $result = $schema->resultset('Seq')->find($gi);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->warn("no result!");
	}
	return $result;
}

=item search_seq

Given a search clause, returns the matching L<Bio::Phylo::PhyLoTA::DAO::Result::Seq> 
objects.

=cut

sub search_seq {
	my ( $self, $clause ) = @_;
	my $result;
	eval {
		$result = $schema->resultset('Seq')->search($clause);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->warn("no result!");
	}
	return $result;	
}

=item single_seq

Given a search clause, returns the matching single 
L<Bio::Phylo::PhyLoTA::DAO::Result::Seq> object.

=cut

sub single_seq {
	my ( $self, $clause ) = @_;
	my $result;
	eval {
		$result = $schema->resultset('Seq')->single($clause);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->warn("no result!");
	}
	return $result;	
}

=item find_node

Given a taxon ID, returns the matching L<Bio::Phylo::PhyLoTA::DAO::Result::Node> object.

=cut

sub find_node {
	my ( $self, $ti ) = @_;
	my $result;
	eval {
		$result = $schema->resultset('Node')->find($ti);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->info("no result!");
	}
	return $result;	
}

=item search_node

Given a search clause, returns the matching L<Bio::Phylo::PhyLoTA::DAO::Result::Node> 
objects.

=cut

sub search_node {
	my ( $self, $clause ) = @_;
	my $result;
	eval {
		$result = $schema->resultset('Node')->search($clause);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->info("no result!");
	}
	return $result;	
}

=item single_node

Given a search clause, returns the single matching 
L<Bio::Phylo::PhyLoTA::DAO::Result::Node> object.

=cut

sub single_node {
	my ( $self, $clause ) = @_;
	my $result;
	eval {
		$result = $schema->resultset('Node')->single($clause);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->info("no result!");
	}
	return $result;	
}

=item search_feature

Given a search clause, returns the matching 
L<Bio::Phylo::PhyLoTA::DAO::Result::Feature> objects.

=cut

sub search_feature {
	my ( $self, $clause ) = @_;
	my $result;
	eval {
		$result = $schema->resultset('Feature')->search($clause);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->info("no result!");
	}
	return $result;
}

=item single_feature

Given a search clause, returns the single matching 
L<Bio::Phylo::PhyLoTA::DAO::Result::Feature> object.

=cut

sub single_feature {
        my ( $self, $clause ) = @_;
        my $result;
        eval {
              	$result = $schema->resultset('Feature')->single($clause);
        };
	if ( $@ ) {
                throw 'BadArgs' => $@;
        }
	if ( not $result ) {
                $logger->info("no result!");
        }
	return $result;
}

=item search_cluster

Given a search clause, returns the matching 
L<Bio::Phylo::PhyLoTA::DAO::Result::Cluster> objects.

=cut

sub search_cluster {
	my ( $self, $clause ) = @_;
	my $result;
	eval {
		$result = $schema->resultset('Cluster')->search($clause);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->info("no result!");
	}
	return $result;	
}

=item single_cluster

Given a search clause, returns the single matching 
L<Bio::Phylo::PhyLoTA::DAO::Result::Cluster> object.

=cut

sub single_cluster {
	my ( $self, $clause ) = @_;
	
	# this hack is necessary because of the strange way the cluster and ci_gi
	# tables are designed. it would have been so much better if there was a
	# simple primary/foreign key relationship, but instead there is a compound
	# key of cluster.ti_root,cluster.ci,cluster.cl_type that matches
	# cigi.ti,cigi.clustid,cigi.cl_type
	if ( exists $clause->{'clustid'} ) {
		my $value = $clause->{'clustid'};
		$clause->{'ci'} = $value;
		delete $clause->{'clustid'};
		$logger->info("search clause included 'clustid', changed this to 'ci'");
	}
	if ( exists $clause->{'ti'} ) {
		my $value = $clause->{'ti'};
		$clause->{'ti_root'} = $value;
		delete $clause->{'ti'};
		$logger->info("search clause included 'ti', changed this to 'ti_root'");
	}	
	
	my $result;
	eval {
		$result = $schema->resultset('Cluster')->single($clause);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->warn("no result!");
	}
	return $result;	
}

=item search_ci_gi

Given a search clause, returns the matching 
L<Bio::Phylo::PhyLoTA::DAO::Result::CiGi> objects.

=cut

sub search_ci_gi {
	my ( $self, $clause ) = @_;
	
	# this hack is necessary because of the strange way the cluster and ci_gi
	# tables are designed. it would have been so much better if there was a
	# simple primary/foreign key relationship, but instead there is a compound
	# key of cluster.ti_root,cluster.ci,cluster.cl_type that matches
	# cigi.ti,cigi.clustid,cigi.cl_type
	if ( exists $clause->{'ci'} ) {
		my $value = $clause->{'ci'};
		$clause->{'clustid'} = $value;
		delete $clause->{'ci'};
		$logger->info("search clause included 'ci', changed this to 'clustid'");
	}
	if ( exists $clause->{'ti_root'} ) {
		my $value = $clause->{'ti_root'};
		$clause->{'ti'} = $value;
		delete $clause->{'ti_root'};
		$logger->info("search clause included 'ti_root', changed this to 'ti'");
	}
	
	my $result;
	eval {
		$result = $schema->resultset('CiGi')->search($clause);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->warn("no result!");
	}
	return $result;	
}

=item search_inparanoid

Given a search clause, returns the matching 
L<Bio::Phylo::PhyLoTA::DAO::Result::Inparanoid> objects.

=cut

sub search_inparanoid {
	my ( $self, $clause ) = @_;
	my $result;
	eval {
		$result = $schema->resultset('Inparanoid')->search($clause);
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->warn("no result!");
	}
	return $result;	
}

=item max_ti

Returns the highest taxon ID from the nodes table.

=cut

sub max_ti {
	my $self = $_;
	my $result;
	eval { 
		$result = $schema->resultset('Node')->get_column('ti')->max;
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->warn("no result!");
	}
	return $result;		
}

=item max_gi

Returns the highest GI from the seqs table.

=cut

sub max_gi {
	my $self = $_;
	my $result;
	eval { 
		$result = $schema->resultset('Seq')->get_column('gi')->max;
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->warn("no result!");
	}
	return $result;		
}

=item insert_seq

Inserts a sequence in the database.

=cut

sub insert_seq {
	my ( $self, $clause ) = @_;
	my $result;
	my @cols = keys(%$clause);
	use Data::Dumper;
	print Dumper(\@cols);
	eval { 
		$result = $schema->populate('Seq', [ \@cols, [@{$clause}{@cols}] ] );
	};
	if ( $@ ) {
		throw 'BadArgs' => $@;
	}
	if ( not $result ) {
		$logger->warn("no result!");
	}
	return $result;			
}

=back

=cut

1;
