# this is an object oriented perl module
package Bio::Phylo::PhyLoTA::DBH;
use strict;
use warnings;
use DBI;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::Logger ':levels';

my $SINGLETON;
our $AUTOLOAD;
my $log = Bio::Phylo::Util::Logger->new;

=head1 NAME

Bio::Phylo::PhyLoTA::DBH - persistent singleton database handle

=head1 SYNOPSIS

 use Bio::Phylo::PhyLoTA::DBH;
 my $dbh = Bio::Phylo::PhyLoTA::DBH->new;
 my $result = $dbh->do('select foo from bar;');

=head1 DESCRIPTION

This package provides a database handle that configures itself based on the values
provided by L<Bio::Phylo::PhyLoTA::Config>, such as database name, user and host. 
Note that no password is used for database connection, since supersmart runs on
a virtual host using a local copy of the public phylota database, therefore privacy
is not an issue. The handle tries to stay-alive, and reconnects to the database if the 
connection goes down. The handle is a singleton object, the idea being that this 
might prevent concurrency issues when running parallel jobs (whether this is the 
right approach is not determined). The handle is used by the L<Bio::Phylo::PhyLoTA::DAO> 
package, there is no direct usage of it elsewhere in the code. 

=head1 METHODS

=over

=item new

The constructor takes the following optional, named arguments:

 -rdbms    => <database system, e.g. mysql> 
 -database => <database name>
 -host     => <location of database server, e.g. localhost>
 -user     => <user name>
 -dsn      => <dsn string>
 -dbh      => <underlying handle>

=back

=cut

sub new {
    my $class = shift;
    my %args  = @_;
    if ( not $SINGLETON ) {
        my $config   = Bio::Phylo::PhyLoTA::Config->new;
        $args{'-rdbms'}    ||= $config->RDBMS;
        $args{'-database'} ||= $config->DATABASE;
        $args{'-host'}     ||= $config->HOST;
        $args{'-user'}     ||= $config->USER;
	my $dsn_tmpl  = 'DBI:%s:database=%s;host=%s';
        $args{'-dsn'} = sprintf($dsn_tmpl, @args{qw[-rdbms -database -host]});
        $args{'-dbh'} = DBI->connect($args{'-dsn'},@args{qw[-user]}, undef, { RaiseError => 1 });
        $SINGLETON = \%args;
        bless $SINGLETON, $class;
    }
    return $SINGLETON;
}

sub DESTROY {
    my $self = shift;
    $self->disconnect if $self->dbh;
}

sub AUTOLOAD {
    my $self = shift;
    my $method = $AUTOLOAD;
    $method =~ s/.+://;
    if ( exists $self->{"-$method"} ) {
        return $self->{"-$method"};
    }
    else {
        if ( not $self->dbh->ping ) {
            $log->warn("handle was disconnected, reconnecting...");
            $self->{'-dbh'} = DBI->connect($self->dsn,$self->user, undef, { RaiseError => 1 });
        }
        $self->dbh->$method(@_);
    }    
}

1;
