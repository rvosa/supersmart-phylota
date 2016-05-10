package Bio::SUPERSMART::App::smrtutils;

use Bio::SUPERSMART::Config;
my $config 	= Bio::SUPERSMART::Config->new;
our $VERSION = $config->RELEASE;

use App::Cmd::Setup -app;


=head1 NAME

Bio::SUPERSMART::App::smrt app class

=head1 DESCRIPTION

This class is invoked when using the command C<smrt>.

=cut

1;
