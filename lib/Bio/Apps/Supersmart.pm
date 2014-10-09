package Bio::Apps::Supersmart;

use App::Cmd::Setup -app;


=head1 NAME

Bio::Apps::Supersmart - 

=head1 DESCRIPTION

=cut

1;

package Bio::Apps::GlobalCmd;

=head1 METHODS

=over

=item opt_spec

Speciefies the global options for the smrt command. This subroutine will
also call the 'options' for the child classes to capture also the
subcommand specific options.

=cut 

sub opt_spec {
	my ($class, $app) = @_;		
	return (
		["verbose|v+", "increase verbosity level", {}],
		["workdir|w=s", "rank to which root taxa are expanded", { default => ".", arg => "dir"}],	
		$class->options($app),
	);	
}

=item usage_desc

Overrides the method from App:Cmd to compile a custom 'usage' 
string for a given subcommand with options and arguments 
of the form 'smrt subcommand option1 <arg1> [option2 <arg2>] [option3]'.  
Options that require arguments followed by the argument name enclosed in '<>'.
Options that are not mandatory are enclosed in []. 

=cut 

sub usage_desc {
	my ($class, $opt, $args) = @_;	
	my $cmd = $class->command_names;	
	my @opts = $class->opt_spec;			
	my $usage = "%c $cmd ";
		
	# build custom usage string (default was "%c $cmd %o" )
	foreach my $opt (@opts){
		my $s = @{$opt}[0];
		my %h = %{@{$opt}[2]};
		my ($short_opt) = ($s =~ /\|([a-z]+)/);
		my $required = $h{'required'};
		my $arg = $h{'arg'};
		my $opt_str = "-$short_opt ";
		$opt_str .= $arg ? " <$arg>"  :"";
		$opt_str = !$required ? "[$opt_str] " : " $opt_str ";		
		$usage .= $opt_str;		
	}
	return $usage;	
}

1;
