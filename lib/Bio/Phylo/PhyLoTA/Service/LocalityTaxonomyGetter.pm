# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Service::LocalityTaxonomyGetter;
use strict;
use warnings;
use Moose;
use IO::File;
use XML::Twig;
use Data::Dumper;
use Bio::Phylo::PhyLoTA::Service::GbifReader;

extends 'Bio::Phylo::PhyLoTA::Service';

=head1 TITLE

Bio::Phylo::PhyLoTA::Service::LocalityTaxonomyGetter - getter of taxon
occurrence data

=head1 DESCRIPTION

This module manages the logic for importing occurrence data from GBIF. This
needs to be done in two steps: i) align the NCBI taxonomy with GBIF's
representation of it. This is done by downloading the spreadsheet that has
the mapping from here: http://data.gbif.org/datasets/resource/13565 and
processing it with store_gbif_ids. ii) fetch the occurrence records from
GBIF's RESTful web service. This is done by passing a filename or URL into
store_gbif_occurrences such as
http://data.gbif.org/ws/rest/occurrence/list?taxonconceptkey=106015799

XXX: this package will be retired.

=head1 METHODS

=over

=cut

# URL stem for requests of occurrence records
my $GBIF_REST = 'http://data.gbif.org/ws/rest/occurrence/list?taxonconceptkey=';

=item store_gbif_occurrences

XXX: irrelevant.

=cut

sub store_gbif_occurrences {
    my ( $self, $resource ) = @_;
    
    # instantiate XML parser
    my $twig = XML::Twig->new(
        'twig_handlers' => {
            
            # set a handler for each occurrence record
            'to:TaxonOccurrence' => sub {
                my ( $twig, $elt ) = @_;
                
                # fetch the GBIF taxon identifier, xpath relative to focal elt
                my $tc = $elt->
                    first_child('to:identifiedTo')->
                    first_child('to:Identification')->
                    first_child('to:taxon')->
                    first_child('tc:TaxonConcept');
                
                # create record representation
                my %record = (
                    'occurrence_id' => $elt->att('gbifKey'),
                    'gbif_id'       => $tc->att('gbifKey'),
                );
                
                # iterate over useful metadata
                for my $key ( qw(country decimalLatitude decimalLongitude) ) {
                    my $child  = $elt->first_child("to:$key");
                    my $value  = $child->text if $child;
                    my $column = $key;
                    $column =~ s/decimalL/l/;
                    $record{$column} = $value;
                }
                
                # let's skip over "empty" records
                return unless $record{'country'} or $record{'latitude'} && $record{'longitude'};                
                $self->logger->debug(Dumper(\%record));
            }
        }
    );
    
    # handle either file or url
    my $method = ( $resource =~ /^http:/ ) ? 'parseurl' : 'parsefile';
    $self->logger->info("going to read from $resource");
    $twig->$method($resource);
    $self->logger->info("done reading from $resource");
}

=item store_gbif_ids

XXX: irrelevant

=cut

sub store_gbif_ids {
    my ( $self, $file, $fetch_occurrences ) = @_;
    
    # instantiate reader 
    my $r = Bio::Phylo::PhyLoTA::Service::GbifReader->new( 'file' => $file );    
    
    # iterate over GBIF species, try to match with NCBI
    while ( my $species = $r->next_species ) {
        my $gbif_id = $species->id;
        
        # search on name initially
        if ( my $gbif_name = $species->binomial ) {
            my @nodes = $self->schema->resultset('Node')->search( { 'taxon_name' => $gbif_name } )->all;
            
            # create a lookup table for names in the GBIF lineage,
            # with their 'depth' (node distance from tip)
            my $depth = 0;
            my %classification = map { $_ => $depth++ } $species->classification;
            
            # iterate over NCBI nodes that match the binomial
            MATCH: for my $node ( @nodes ) {
                my $ti = $node->ti;
                
                # traverse classification for focal NCBI node
                my ( $seen, %lookup ) = ( 0 );
                NODE: while( $node ) {
                    last NODE if $seen == $depth;
                    my $name = $node->taxon_name;
                    
                    # these don't all exist because NCBI has more
                    # levels than GBIF
                    if ( exists $classification{$name} ) {
                        
                        # if $seen is higher, it means levels were
                        # in different order
                        if ( $seen <= $classification{$name} ) {
                            $lookup{$name} = $seen++;
                        }
                        else {
                            $self->logger->warn("apparent mismatch in higher ranks of $gbif_name ($ti)");
                            next MATCH;
                        }
                    }
                    
                    # go to parent
                    $node = $self->schema->resultset('Node')->find($node->ti_anc);
                }
                
                # if we reach this point, we've found a valid match
                $self->logger->info("GBIF:$gbif_id <=> NCBI:$ti ($gbif_name)");
                $self->logger->debug("GBIF:".Dumper(\%classification));
                $self->logger->debug("NCBI:".Dumper(\%lookup));
                
                # optionally follow links to occurrence records
                if ( $fetch_occurrences ) {
                    $self->store_gbif_occurrences( $GBIF_REST . $gbif_id );
                }
            }
        }
        else {
            $self->logger->warn("no binomial for GBIF record " . $species->id());
        }
    }
}

=back

=cut

1;
