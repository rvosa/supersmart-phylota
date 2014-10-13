package Bio::Apps::Supersmart;
our $VERSION = 0.1;
use App::Cmd::Setup -app;

=head1 NAME

Bio::Apps::Supersmart app class

=head1 DESCRIPTION

This class is invoked when using the command C<smrt>.

=cut

1;

package Bio::Apps::GlobalCmd;

=head1 NAME

Bio::Apps::GlobalCmd Superclass for C<smrt> subcommands

=head1 DESCRIPTION

This class provides global options (e.g. verbosity, working directory) that are
inherited by all subcommand classes.

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
		[ "help|h", "display help screen", {} ],
		[ "verbose|v+", "increase verbosity level", {} ],
		[ "workdir|w=s", "directory in which results and intermediate files are stored", { default => ".", arg => "dir"} ],	
		$class->options($app),
	);	
}

=item description

Returns the description for a specific subcommand. This method is invoked by 
the child classes when the global option 'help' is used, e.g. in 
$ smrt help subcommand
The actual description of the subcommand is parsed from it's DESCRIPTION field
in the POD documentation.

=cut

sub description {
	my ($class, $opt, $args) = @_;

	# classname to filename
  	my $pm_file = $class;
  	$pm_file =~ s/::/\//g; 
  	$pm_file =~ s/=.*//g;   	
  	$pm_file .= '.pm';
  	$pm_file = $INC{$pm_file} or return "(unknown)";
	
  	# if the pm file exists, open it and parse it
  	open my $fh, "<", $pm_file or return "(unknown)";
	undef $/;	
	my $content = <$fh>;
	my ($desc) = 
         ($content =~ m/.*?DESCRIPTION\n(.*?)\n=.*/s);	
	close $fh;
	$content =~ s/^\s+//;	
	return $desc;
}

=item validate_args

overrides the default method mainly for one reason: There are two ways of 
displaying the help screen: One is '$ smrt help command', the other is
'$ smrt command -h'. Here we make the latter to behave exactly as
the former.

=cut

sub validate_args {
	my ($class, $opt, $args) = @_;		
	
	if ($opt->help){
		# This is a small hack to make the help screen called with the 
		#  option (as here, e.g. "$ smrt command -h" ) look the same
		#  as when we call "$ smrt help command"  				
		my ($cmd) = $class->app->prepare_command("help"); 
		my @name = ($class->command_names);
		$cmd->execute($opt, \@name);
		exit;
	} else {
		$class->validate($opt, $args);
	}
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
		my $required = $h{'mandatory'};
		my $arg = $h{'arg'};
		my $opt_str = "-$short_opt ";
		$opt_str .= $arg ? " <$arg>"  :"";
		$opt_str = !$required ? "[$opt_str] " : " $opt_str ";		
		$usage .= $opt_str;		
	}
	return $usage;	
}

=back

1;
