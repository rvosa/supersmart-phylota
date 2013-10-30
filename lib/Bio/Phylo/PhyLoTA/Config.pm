# this is an object oriented perl module
package Bio::Phylo::PhyLoTA::Config;
use strict;
use warnings;
use Config::Tiny;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(printval);
our $AUTOLOAD;
my $SINGLETON;

sub new {
    if ( not $SINGLETON ) {
        return $SINGLETON = shift->read(@_);    
    }
    else {
        return $SINGLETON;
    }
}

sub printval {
	my $key  = shift @ARGV;
	my $self = __PACKAGE__->new;
	print $self->$key;
}

sub read {
    my $self = shift;
    my $file = shift || "$ENV{PHYLOTA_HOME}/conf/phylota.ini";
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

sub currentConfigFile { shift->{'_file'} }

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
        
        # make paths absolute if that resolves an otherwise non-existant path
        if ( $key =~ /(?:FILE|DIR)$/ ) {
            if ( not -e $root->{$key} and -e $ENV{PHYLOTA_HOME} . '/' . $root->{$key} ) {
                return $ENV{PHYLOTA_HOME} . '/' . $root->{$key};
            }            
        }        
        return $root->{$key};
    }
    else {
        warn "Unknown key: $key" unless $key =~ /^[A-Z]+$/;
    }
}

1;