# this is an object oriented perl module
package Bio::Phylo::PhyLoTA::Config;
use strict;
use warnings;
use Config::Tiny;
use Bio::SUPERSMART;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(printval);
our $AUTOLOAD;
our $VERSION = $Bio::SUPERSMART::VERSION;
my $SINGLETON;

=head1 NAME

Bio::Phylo::PhyLoTA::Config - manager of runtime configuration variables

=head1 SYNOPSYS

 use Bio::Phylo::PhyLoTA::Config; 
 
 # instantiate the singleton config object
 my $conf = Bio::Phylo::PhyLoTA::Config->new;
 
 # fetch the value of a variable in supersmart.ini
 my $examl = $conf->EXAML_BIN;

=head1 DESCRIPTION

The config object holds the values of variables that are needed when running the
pipeline. Most commonly, these are the locations of executables that the pipeline
invokes. The values of these variables are defined in C<conf/supersmart.ini>, although
you can override these using environment variables. For example, to define an
alternate location for EXAML_BIN, create an environment variable called 
SUPERSMART_EXAML_BIN, whose value will be used in lieu of that in the config file.

=head1 METHODS

=over

=item new

The constructor optionally takes the location of a configuration file in INI
syntax. By default this is C<conf/supersmart.ini>.

=cut

sub new {
    if ( not $SINGLETON ) {
        $SINGLETON = shift->read(@_);
		
		# To allow for reproducibility of stochastic steps, such as
		# bootstrapping, random topology resolution, tree generation, markov
		# chains, and so on, we have defined a random seed. This gets passed
		# to the various programs that allow for such a seed to be provided.
		# In addition to that, we need to seed Perl's random number generator
		# for the instances where we do stochastic things in our own code.
		# As this seeding must happen only once, doing it here inside the
		# pseudo-constructor of this singleton object is a decent place.
		if ( $SINGLETON->RANDOM_SEED ) {
			srand $SINGLETON->RANDOM_SEED;
		}
		return $SINGLETON;
    }
    else {
        return $SINGLETON;
    }
}

=item printval

Prints the value of a variable to STDOUT. Example on the shell:

 $ perl -MBio::Phylo::PhyLoTA::Config=printval -e 'printval "EXAML_BIN"'

=cut

sub printval {
	my $key  = shift @ARGV;
	my $self = __PACKAGE__->new;
	print $self->$key;
}

=item read

Reads a configuration file in INI syntax, populates the config object 
with the values therein. Optional argument is an alternate location
for the config file.

=cut

sub read {
    my $self = shift;
    my $file = shift || "$ENV{SUPERSMART_HOME}/conf/supersmart.ini";
    my $conf = Config::Tiny->read($file);
    if ( my $error = Config::Tiny->errstr ) {
        die $error;
    }
    if ( ref $self ) {
        $self->{'_conf'} = $conf;
    }
    else {
        my $class = $self;
        $self = bless { '_conf' => $conf }, $class;
    }
    $self->{'_file'} = $file;
    return $self;
}

=item currentConfigFile

Returns the location of the currently used configuration file.

=cut

sub currentConfigFile { shift->{'_file'} }

=item currentGBRelease

Returns the currently used GenBank release number (i.e., 184).

=cut

sub currentGBRelease {
    my $self = shift;
    if ( ! $self->GB_RELNUM ) {
        my $file = $self->GB_RELNUM_FILE;
        open my $fh, '<', $file or die $!;
        while(<$fh>) {
            chomp;
            $self->GB_RELNUM($_) if $_;
        }
        close $fh;
    }
    return $self->GB_RELNUM;
}

=item currentGBReleaseDate

Returns the date stamp for the currently used GenBank release.

=cut

sub currentGBReleaseDate {    
    my $self = shift;
    if ( ! $self->GB_RELNUM_DATE ) {
        my $file = $self->GB_RELNUM_DATE_FILE;
        open my $fh, '<', $file or die $!;
        while(<$fh>) {
            chomp;
            $self->GB_RELNUM_DATE($_) if $_;
        }
        close $fh;
    }
    return $self->GB_RELNUM_DATE;
}

=item RELEASE

Returns the release number.

=cut

sub RELEASE { $VERSION }

sub AUTOLOAD {
    my $self = shift;
    my $root = $self->{'_conf'}->{'_'};
    my $key = $AUTOLOAD;
    $key =~ s/.+://;
    if ( exists $ENV{"SUPERSMART_${key}"} ) {
    	return $ENV{"SUPERSMART_${key}"};
    }
    if ( exists $root->{$key} ) {
        $root->{$key} = shift if @_;
		
		# dereference environment variables if present        
        $root->{$key} =~ s/\$\{(\w+)\}/$ENV{$1}/g;        
        
        # make paths absolute if that resolves an otherwise non-existant path
        if ( $key =~ /(?:FILE|DIR)$/ ) {
            if ( not -e $root->{$key} and -e $ENV{SUPERSMART_HOME} . '/' . $root->{$key} ) {
                return $ENV{SUPERSMART_HOME} . '/' . $root->{$key};
            }            
        }        
        return $root->{$key};
    }
    else {
        warn "Unknown key: $key" unless $key =~ /^[A-Z]+$/;
    }
}

=back

=cut

1;
