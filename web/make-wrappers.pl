use strict;
use warnings;

use Data::Dumper;
use Module::Find;

use Bio::SUPERSMART::Config;

# get all subcommands that have a subcommand
my @modules = usesub(Bio::SUPERSMART::App::smrt::Command);
my @subcommands = map { $_=~ s/.*:://g; $_} @modules;

my $subcommand = "Taxize";

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

# given a subcommand, constructs the galaxy 'tool' tag and returns it as a hash
sub get_tool {
	my $subcommand = shift;
	
	my %result;

    # gather attributes: id, name and version
	my $c = Bio::SUPERSMART::Config->new;	
	my $version = $c->RELEASE;	
	$result{'tool'} = {'id' => "smrt_$subcommand", 'name' => $subcommand, 'version' => "$version"};

	# get and set description
	my $description = "Bio::SUPERSMART::App::smrt::Command::$subcommand"->description;
	$result{'tool'}->{'description'} = [ $description ];

	# prevent stderr from being error in galaxy:
	$result{'tool'}->{'stdio'}->{'exit_code'} = { 'range' => '1:', 'err_level' => 'fatal'};
	
	my $in_out = get_inputs_outputs($subcommand);
	my @a = values(%{$in_out});

	@result{keys(%{$in_out})} = values(%{$in_out});

	print Dumper(\%result);

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
	my @in = grep {$_->{'param'}} @in_out;
	my @out = grep {$_->{'data'}} @in_out;
	
	# put into higher level structure
	$result{'inputs'} = \@in if scalar @in;
	$result{'outputs'} = \@out if scalar @out;
	
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

	# only process when option is desired to appear in Galaxy,
	# as set by the attributes galaxy_in and galaxy_out in the command class
	if ( $info{'galaxy_in'} || $info{'galaxy_out'} ) {
		$h{'name'} = $name_str;
		$h{'label'} = $name_str;		
	} 
	else {
		return ();
	}
	   	
	my %result = ( $tag => \%h );
	return \%result;
}



