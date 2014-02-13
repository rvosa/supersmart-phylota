package Bio::Phylo::PhyLoTA::Service::GbifReader;
use Moose;
use Bio::Species;
use IO::File;

has 'header' => (
	isa => 'ArrayRef[Str]',
	is  => 'rw'	
);

has 'handle' => (
	isa => 'IO::Handle',
	is  => 'rw',
	'trigger' => sub {
		my ( $self, $handle ) = @_;
		if ( $handle ) {
			my @header = split /\t/, $handle->getline;
			$self->header(\@header);
		}
	}
);

has 'file' => (
	isa => 'Str',
	is  => 'rw',
	'trigger' => sub {
		my ( $self, $file ) = @_;
		if ( $file ) {
			my $fh = IO::File->new;
			$fh->open( $file, '<' );
			$self->handle($fh);
		}
	}
);

=head1 NAME

Bio::Phylo::PhyLoTA::Service::GbifReader - reads GBIF data

=head1 DESCRIPTION

XXX: this package will be retired.

=head1 METHODS

=over

=item next_species

Returns the next L<Bio::Species> object from the underlying GBIF data file.

=back

=cut

sub next_species {
	my $self = shift;
	
	# attempt to read the next line
	my @fields = split /\t/, $self->handle->getline;
	
	# we may be at the end of the file, so this could be false
	if ( @fields ) {
		
		# create a hash with keys from the file header
		my @header = @{ $self->header };
		my %record = map { $header[$_] => $fields[$_] } 0 .. $#header;
		
		# order the fields as Bio::Species' @classification array expects
		my @classfields = qw(Species Genus Family Order Class Phylum Kingdom);
		
		# get the values
		my %seen;
		my @class = grep { ! $seen{$_}++ } grep { /\S/ } @record{@classfields};
		
		# done
		return Bio::Species->new(
			'-classification' => \@class,
			'-species'        => $record{'Species'},
			'-id'             => $record{'GBIF ID'},
		);
	}
	return undef;
}

1;