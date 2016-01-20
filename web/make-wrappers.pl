use strict;
use warnings;

use Data::Dumper;
use Module::Find;
use XML::Simple;

use Bio::SUPERSMART::Config;

# get all subcommands that have a subcommand
my @modules = usesub(Bio::SUPERSMART::App::smrt::Command);
my @subcommands = map { $_=~ s/.*:://g; $_} @modules;

my $subcommand = "Align";

#get_inputs_outputs($subcommand);
get_tool($subcommand);




sub to_galaxy_xml {
	
	# tool tag
	
	# description tag

	# command

	# inputs

	# outputs

	# help

}

# given a subcommand, constructs the galaxy 'tool' tag (the root tag for a tool) and returns it as a hash
sub get_tool {
	my $subcommand = shift;
	
	# hash representing xml structure of entire tool
	my %result;
	
	# string containing detailed information, will be in <help> tag
	my $help;

    # gather attributes: id, name and version
	my $c = Bio::SUPERSMART::Config->new;	
	my $version = $c->RELEASE;	
	$result{'tool'} = {'id' => "smrt_$subcommand", 'name' => $subcommand, 'version' => "$version"};

	# get and set description
	my $description = "Bio::SUPERSMART::App::smrt::Command::$subcommand"->description;
	my $abstract = "Bio::SUPERSMART::App::smrt::Command::$subcommand"->abstract;
	$description =~ s/\s+$//g;	
	$help .= "**What it does**\n\n$description\n";	
	$result{'tool'}->{'description'} = [ "\n$abstract\n" ];
	

	# prevent stderr from being error in galaxy:
	$result{'tool'}->{'stdio'}->{'exit_code'} = [ { 'range' => '1:', 'err_level' => 'fatal'} ];
	
	# set in- and output parameters
	my $in_out = get_inputs_outputs($subcommand);

	# add workspace in- and outputs
	$result{'tool'}->{'inputs'} = $in_out->{'inputs'};
	$result{'tool'}->{'outputs'} = $in_out->{'outputs'};

	# make command tag
	my $cmd; 
	
	# unzip workspace, if exists
	#$cmd .= "\n#if \$workspace\n";
	#$cmd .= "\tunzip \$workspace;\n";
	#$cmd .= "#end if\n\n";
	#
	$cmd .= "smrt " . lc($subcommand) . "\n";
	
	# add input and output arguments to command
	my @inputs = @{ $in_out->{'inputs'}->{'param'} };
	my @outputs = @{ $in_out->{'outputs'}->{'data'} };
	for ( @inputs ) {
		my %h = %{$_};
		my $param_name = $h{"name"};
		$cmd .= "#if \$${param_name}\n";
		$cmd .= "\t--$param_name \$${param_name}\n";
		$cmd .= "#end if\n\n";
	}
	for ( @outputs ) {
		my %h = %{$_};
		my $param_name = $h{"name"};

		$cmd .= "\t--$param_name \$${param_name}\n\n" unless $param_name eq "workspace";		
	}

	# add workidir
	$cmd .= "--workdir /tmp/\$jobid;\n\n";

	# zip workspace which is returned by galaxy
	#$cmd .= "zip workspace.zip *; mv workspace.zip \$workspace;\n\n";

	# set command
	$result{'tool'}->{'command'} = [ $cmd ];

	# set help
	$result{'tool'}->{'help'} = [ $help ];

	my $out = XMLout(\%result, KeepRoot => 1);
	print $out;
	
	return \%result;

}



# given a command, gets the tags 'inputs' and 'outputs' for galaxy xml as hashrefs
sub get_inputs_outputs {
	my $subcommand = shift;

	my %result;

	# collect options 
	my @options = "Bio::SUPERSMART::App::smrt::Command::$subcommand"->options;
	
	# get hashs with input and outputs from subcommand
	my @in_out = map { parse_option($_) } @options;
	
	# seperate in- and output elements into a higher level structure,
	# put into tags 'inputs' and 'outputs'. In Galaxy, per default 
	# inputs have the tag 'param' and outputs the tag 'data'
	my @in = map { values %$_ } grep { $_->{'param'} } @in_out;
	my @out = map { values %$_ } grep { $_->{'data'} } @in_out;

	my @kk = map {values %$_} @in;
	
	# put into higher level structure
	$result{'inputs'}->{'param'} = \@in if scalar @in;
	$result{'outputs'}->{'data'} = \@out if scalar @out;
	
	# add workspace for input data
	push $result{'inputs'}->{'param'}, {'name'=>'jobid', 'label'=>'jobid', 'type'=>'text'};
	#$result{'inputs'}->{'data'} = [ {'name'=>'workspace', 'label'=>'workspace', 'type'=>'zip'} ];
	#push $result{'outputs'}->{'data'}, {'name'=>'workspace', 'label'=>'workspace', 'type'=>'zip'};

	return \%result;
}

# parses single App::Cmd command-line option and returns hash with options
sub parse_option {
	my $op = shift;

	my @arr = @{$op};

	my $name_str = $arr[0];
	my $description = $arr[1];

	my %info = %{$arr[2]};
	
	# extract long option name  from name string
	$name_str =~ s/\|.+$//g;
	
	# set tag for xml: per default in Galaxy, for inputs, tag is 'param',
	#  for outputs, tag is 'data'
	my $tag = $info{'galaxy_in'} ? 'param' : 'data';
	my %h;

	# extract type
	my $type = $info{'galaxy_type'};
	if ( not $type ) {
		warn("No galaxy type for option $name_str given. Skipping");
		return ();
	}
	$h{'type'} = $type;

	# extract format
	my $format = $info{'galaxy_format'};
	$h{'format'} = $format if $format;

	# only process when option is desired to appear in Galaxy,
	# as set by the attributes galaxy_in and galaxy_out in the command class
	if ( $info{'galaxy_in'} || $info{'galaxy_out'} ) {
		$h{'name'} = $name_str;
		$h{'label'} = "$name_str: $description";		
	} 
	else {
		return ();
	}
	   	
	my %result = ( $tag => \%h );
	return \%result;
}

#parameter types:
=pod
   text=TextToolParameter,
    integer=IntegerToolParameter,
    float=FloatToolParameter,
    boolean=BooleanToolParameter,
    genomebuild=GenomeBuildParameter,
    select=SelectToolParameter,
    color=ColorToolParameter,
    data_column=ColumnListParameter,
    hidden=HiddenToolParameter,
    hidden_data=HiddenDataToolParameter,
    baseurl=BaseURLToolParameter,
    file=FileToolParameter,
    ftpfile=FTPFileToolParameter,
    data=DataToolParameter,
    data_collection=DataCollectionToolParameter,
    library_data=LibraryDatasetToolParameter,
    drill_down=DrillDownSelectToolParameter
=cut
