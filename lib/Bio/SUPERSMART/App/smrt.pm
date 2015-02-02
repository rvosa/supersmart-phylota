package Bio::SUPERSMART::App::smrt;

use Bio::Phylo::PhyLoTA::Config;
my $config 	= Bio::Phylo::PhyLoTA::Config->new;
our $VERSION = $config->RELEASE;

use App::Cmd::Setup -app;


=head1 NAME

Bio::SUPERSMART::App::smrt app class

=head1 DESCRIPTION

This class is invoked when using the command C<smrt>.

=cut

1;
