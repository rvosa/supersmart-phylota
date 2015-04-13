package Bio::SUPERSMART::App::smrt::SubCommand;

use Cwd;
use Bio::Phylo::Util::Logger ':levels';

=head1 NAME

Bio::SUPERSMART::App::smrt::SubCommand Superclass for C<smrt> subcommands

=head1 DESCRIPTION

This class provides global options (e.g. verbosity, working directory) that are
inherited by all subcommand classes.

=cut

my $logfh;

=head1 METHODS

All methods below are inherited by child classes.

=over

=item logger

Getter/setter for the subcommand's Bio::Phylo::Util::Logger object, 
verbosity can be de- or increased by  the argument -v for all subcommands. 

=cut

sub logger {
	my $self = shift;
	if ( @_ ) {
		$self->{'logger'} = shift;		
	}
	return $self->{'logger'};
}

=item outfile

returns the output file name of the subcommand

=cut

sub outfile {
	my $self = shift;
	if ( @_ ) {
		$self->{'outfile'} = $self->absolute_path(shift);		
	}
	return $self->{'outfile'};
}

=item workdir

Getter/setter for the working directory for the subcommand

=cut

sub workdir {
	my $self = shift;
	if ( @_ ) {
		$self->{'workdir'} = shift;		
	}
	return $self->{'workdir'};
}

=item absolute_path

creates an absolute path to a file location, if file is not yet
given as absolute path. This is done by prepending the working directory
to the file name.

=cut

sub absolute_path {
	my $self = shift;
	my $filename = shift;
	if (not $filename =~ /\//) {
		if ( my $wd = $self->workdir ) {	
			$filename = $filename =~ /^\// ? $filename : $wd . "/" . $filename;
		} else {
			$self->logger->warn("no working directory specified, using relative paths ")
		}
	}
	return $filename; 
}

=item run

This is the method in which the functionality of the subcommand
class is implemented; all subcommand-specific magig will happen in
here. 

=cut

sub run {	
	die "Subroutine 'run' not implemented for subcommand class " . ref ($_); 		
}

=item execute

Overrides the default classes' 'execute' method such that some global properties 
for all subcommands (e.g. verbosity level, ...) can be set in this method. The
method then calls the 'run' subroutine, which must be implemented for all child classes

=cut

sub execute {
	my ($class, $opt, $args) = @_;
	my $config 	= Bio::Phylo::PhyLoTA::Config->new;
	my $release = $config->RELEASE;
	$class->logger->info("This is SUPERSMART release " . $config->RELEASE);
	
	my $result = $class->run( $opt, $args );	
	close $logfh;
	return $result;
}

=item init

The init subroutine initializes objects that are shared by all subcommands (as for instance the working 
directory, output file and a L<Bio::Phylo::Util::Logger> instance). Also, we check for command line options that require a file argument
and, if necessary, turn their paths into absolute paths.

=cut

sub init {
	my ($self, $opt, $args) = @_;
	my $verbosity = INFO;
	$verbosity += $opt->verbose ? $opt->verbose : 0;
    
	# set working directory
	($wd = $opt->workdir || getcwd()) =~ s/\/$//g;
	$self->workdir($wd);

	# loop through options to see which ones are file options; 
	#  set absolute path for all filenames given
	my %file_opts = map { (my $optname= $_->[0])=~s/\|.+//g; $_->[2]->{'arg'} eq 'file' ? ($optname=>1) : () } $self->options;
	$file_opts{'logfile'} = 1;
	for my $given_opt ( keys %$opt ) {
		if ( exists $file_opts{$given_opt} ) {
			my $new_filename = 
			$opt->{$given_opt} = $self->absolute_path($opt->{$given_opt});
		}		
	}
 	# set outfile name
       if ( my $of = eval { $opt->outfile } ) {
        	$self->outfile($of);
       }
 
 	# create logger object with user-defined verbosity
	$self->logger( Bio::Phylo::Util::Logger->new(
		'-level' => $verbosity,
		'-style' => $opt->logstyle,
		'-class' => [ ref( $self ), 'Bio::SUPERSMART::App::smrt::SubCommand', 
									'Bio::Phylo::PhyLoTA::Service::ParallelService', 
									'Bio::Phylo::PhyLoTA::Service::TreeService', 
									'Bio::Phylo::PhyLoTA::Service::CalibrationService' ],		
    ));
    
    # redirect logger output to file if specified by command call
    if ( my $logfile = $opt->logfile ) {
    	open $logfh, '>', $logfile or die $!;
		$self->logger->set_listeners(sub{$logfh->print(shift)});
    } 
}

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
		[ "workdir|w=s", "directory in which results and intermediate files are stored", { arg => "dir"} ],
		[ "logfile|l=s", "write run-time information to logfile", { arg => "file" }],	
		[ "logstyle|y=s", "toggles logging style between 'simple' and 'detailed'", { default => "simple" }],
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

Overrides the default method mainly for one reason: There are two ways of 
displaying the help screen: One is '$ smrt help command', the other is
'$ smrt command -h'. Here we make the latter to behave exactly as
the former.

=cut

sub validate_args {
	my ($class, $opt, $args) = @_;		
	
	# first init the properties of the superclass
	$class->init($opt, $args);
			
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

=item options

This is a wrapper for L<App::Cmd>'s 'opt_spec' method and is implemente by the child classes to
specify command-line options for a specific subcommand.

=cut

sub options {
	die("subroutine 'options' must be implemented by the respective child class");
}

=item validate

This is a wrapper for L<App::Cmd>'s 'validate_args' method and is implemente by the child classes to
validate command-line options for a specific subcommand.

=cut

sub validate {
	die("subroutine 'validate' must be implemented by the respective child class");

}

=item run

This is a wrapper for L<App::Cmd>'s 'execute' method and is implemente by the child classes to
run a specific subcommand
=cut

sub run {
	die("subroutine 'run' must be implemented by the respective child class");

}

=back

=cut

1;
