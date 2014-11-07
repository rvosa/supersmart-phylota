package Bio::SUPERSMART::App::smrt::SubCommand;

use Bio::Phylo::Util::Logger ':levels';


my $_verbosity = INFO;
my $_logger;
my $_outfile; 
my $_workdir;

=head1 NAME

Bio::SUPERSMART::App::smrt::SubCommand Superclass for C<smrt> subcommands

=head1 DESCRIPTION

This class provides global options (e.g. verbosity, working directory) that are
inherited by all subcommand classes.


=head1 METHODS

All methods below are ingerited by child classes.

=over

=item logger

Getter/setter for the subcommand's Bio::Phylo::Util::Logger object, 
verbosity can be de- or increased by  the argument -v for all subcommands. 

=cut

sub logger {
	my $self = shift;
	if ( @_ ) {
		$_logger = shift;		
	}
	return $_logger;
}

=item outfile

returns the output file name of the subcommand

=cut

sub outfile {
	my $self = shift;
	if ( @_ ) {
		$_outfile = shift;		
	}
	return $_outfile;
}

=item outfile

returns the working directory the subcommand

=cut

sub workdir {
	my $self = shift;
	if ( @_ ) {
		$_workdir = shift;		
	}
	return $_workdir;
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
	return $class->run( $opt, $args );	
}


sub init {
	my ($class, $opt, $args) = @_;
	my $str = ref($class);
	$_verbosity += $opt->verbose ? $opt->verbose : 0;
	
	# create logger object with user-defined verbosity
	$_logger = Bio::Phylo::Util::Logger->new(
		'-level' => $_verbosity,
		'-class' => [ ref( $class ), 'Bio::Phylo::PhyLoTA::Service::ParallelService', 'Bio::Phylo::PhyLoTA::Service::TreeService' ],		
    );
    
    # make output file name. If output file is given and it is an absolute path, 
    #  leave as is; if it is a single filename or a relative path, prepend the working
    #  directory    	
	($wd = $opt->workdir) =~ s/\/$//g;
	$_workdir = $wd;
    my $of = eval { $opt->outfile };
    if ( $of ){
    	$_outfile = $of =~ /^\// ? $of : $wd . "/" . $of;  
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
		[ "workdir|w=s", "directory in which results and intermediate files are stored", { default => "./", arg => "dir"} ],	
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
		
	$class->usage_error("workdir needs to be a valid directory") if not -d $opt->workdir;
	
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

=back

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


1;