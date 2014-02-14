# this is an object oriented perl module
package Bio::Phylo::PhyLoTA::DAO;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

use strict;
use warnings;

use base 'DBIx::Class::Schema';

__PACKAGE__->load_namespaces;


# Created by DBIx::Class::Schema::Loader v0.07002 @ 2012-05-26 14:28:40
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:qy8gwA4ReJ1Jez7iqtYMBQ


# You can replace this text with custom content, and it will be preserved on regeneration

=head1 NAME

Bio::Phylo::PhyLoTA::DAO - The database schema

=head1 SYNOPSIS

    use Bio::Phylo::PhyLoTA::DAO;
    
    my $schema = Bio::Phylo::PhyLoTA::DAO->new;
    my $node = $schema->resultset('Node')->find(9606);
    while($node) {
        print $node->taxon_name, "\n";
        my $clusters = $node->clusters;
        while( my $c = $clusters->next ) {
            my $gis = $schema->resultset('CiGi')->search({ clustid => $c->ci });
            while ( my $g = $gis->next ) {
                print $g->gi, "\n";
            }
        }
        $node = $schema->resultset('Node')->find($node->ti_anc);
    }

=head1 DESCRIPTION

The pipeline uses a relational schema implemented in MySQL with a small number of tables.
These tables are mapped to object classes inside the Bio::Phylo::PhyLoTA::DAO::* 
namespace, principally using the code generator from L<DBIx::Class>. This package
represents the schema as a whole, through which different tables can be queried, as
shown in the SYNOPSIS. Note that this provides a fairly verbose API, and that a more
concise equivalent one is provided by L<Bio::Phylo::PhyLoTA::Service>, which delegates
to this one.

=head1 METHODS

=over

=item new

The constructor returns a singleton object and takes no arguments.

=back

=cut

use Bio::Phylo::PhyLoTA::DBH;
use Bio::Phylo::Util::Logger ':levels';

my $SINGLETON;
my $log = Bio::Phylo::Util::Logger->new;
my %args = ( 'limit_dialect' => 'LimitXY' );

sub new {
	my $package = shift;
	if ( not $SINGLETON ) {
		$log->info("first call to constructor");
	
		# the SUPER::connect method can be passed a code reference
# 		my $sub = sub {
# 			$log->info("executing code ref that returns database handle");
# 			return Bio::Phylo::PhyLoTA::DBH->new;
# 		};		
		my $dbh = Bio::Phylo::PhyLoTA::DBH->new;
		
		# create the singleton
# 		$SINGLETON = $package->connect( $sub, \%args );
		$SINGLETON = $package->connect( $dbh->dsn, $dbh->user, undef, \%args );
	}
	else {
		$log->info("additional, no-op call to singleton constuctor");
	}	
	return $SINGLETON;
}
1;
