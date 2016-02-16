use strict;
use warnings;

use Data::Dumper;
use Module::Find;
use XML::Simple;
use List::MoreUtils qw (uniq);

use Bio::SUPERSMART::Config;

# get all subcommands that have a subcommand
my @modules = usesub(Bio::SUPERSMART::App::smrt::Command);
my @subcommands = map { $_=~ s/.*:://g; $_} @modules;

#@subcommands = ("Taxize", "Align", "Orthologize", "BBmerge"); ##("BBinfer");##reverse ("Taxize", "Align", "Orthologize", "BBmerge", "BBinfer");

for my $subcommand ( @subcommands ) {
	print "Processing $subcommand \n";
	my $tags = get_tool($subcommand);

	my $out = XMLout($tags, KeepRoot => 1, NoEscape => 1);
	my $filename = lc "$subcommand.xml";
	open my $fh, '>',  $filename or die $!;
	print $fh $out;
	close $fh;	
}

##if $condition_inferencetool.inferencetool == "RAxML"
#	--rapid_boot $condition_inferencetool.rapid_boot
##end if


# generate the cheetah code used in the <command> tag
sub _get_cheetah {
	my ($subcommand, $in_out) = @_;
	
	my $cheetah;

	$cheetah .= "smrt " . lc($subcommand) . "\n\n";
	
	# add input parameters to cheetah command
	my @inputs = @{ $in_out->{'inputs'}->{'param'} };
	for ( @inputs ) {
		my %h = %{$_};
		my $param_name = $h{"name"};
		next if $param_name eq 'jobid';
		$cheetah .= "#if \$${param_name}\n";
		$cheetah .= "\t--$param_name \$${param_name}\n";
		$cheetah .= "#end if\n\n";	
	}

	# add output parameters to cheetah command
	if ( $in_out->{'outputs'} ) {
		my @outputs = @{ $in_out->{'outputs'}->{'data'} };
		for ( @outputs ) {
			my %h = %{$_};
			my $param_name = $h{"name"};
			$cheetah .= "--$param_name \$${param_name}\n\n";		
		}
	}

	# add conditional parameters to cheetah command
	if ( $in_out->{'inputs'}->{'conditional'} ) {
		my %conditionals = %{ $in_out->{'inputs'}->{'conditional'} };
		for my $k ( keys %conditionals ) {
			my $c = $conditionals{$k};
			my $condition_name = $c->{'name'};
			
			for my $w ( @{$c->{'when'}} ) {
				my $value = $w->{'value'};
				$cheetah .= "#if \$${condition_name}.${k} == \"${value}\" \n";
				for my $p ( @{$w->{'param'}} ) { 
					my $paramname = $p->{'name'};
					$cheetah .= "\t --$paramname \$${condition_name}.${paramname}\n";
				}
				$cheetah .= "#end if\n\n";
			}		
		}
	}

	# add workidir
	$cheetah .= "--workdir /tmp/\$jobid;\n\n";

	return $cheetah;
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
	my $cmd = _get_cheetah( $subcommand, $in_out ); 
	
	# set command
	$result{'tool'}->{'command'} = [ $cmd ];

	# set help
	$result{'tool'}->{'help'} = [ $help ];

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

	# some parameters have dependencies, they will be wrapped in <conditional> tags
	my $conditionals = get_conditionals(\@in);
	
	# extract all names params that are present in conditionals	and remove them from the input params
	my @ref_params = map { $_->{'param'}->{'name'} } values ( %{$conditionals} );
	my @refd_params = map {$_->{'name'}} map { @{$_->{'param'}} } map { @{$_->{'when'}} } values ( %{$conditionals} );
	my %cond_params = map {$_=>1} uniq ( @ref_params, @refd_params );
	@in = grep { ! $cond_params{$_->{'name'}} } @in;

	# set in, out and conditional parameters 
	$result{'inputs'}->{'conditional'} = $conditionals if keys %{$conditionals};
	$result{'inputs'}->{'param'} = \@in if scalar @in;
	$result{'outputs'}->{'data'} = \@out if scalar @out;
	
	# add workspace for input data
	push $result{'inputs'}->{'param'}, {'name'=>'jobid', 'label'=>'jobid', 'type'=>'text'} if $result{'inputs'}->{'param'};

	return \%result;
}

# parses single App::Cmd command-line option and returns hash with options
sub parse_option {
	my $op = shift;

	my %result;

	my @arr = @{$op};

	my $name_str = $arr[0];
	my $description = $arr[1];
	my %info = %{$arr[2]};

	# only process when option is desired to appear in Galaxy,
	# as set by the attributes galaxy_in and galaxy_out in the command class
	if ( $info{'galaxy_in'} || $info{'galaxy_out'} ) {
	
		# extract long option name from name string
		$name_str =~ s/\|.+$//g;
		
		# set tag for xml: per default in Galaxy, for inputs, tag is 'param',
		#  for outputs, tag is 'data'
		my $tag = $info{'galaxy_in'} ? 'param' : 'data';
		my %h;
		
		# extract type
		my $type = $info{'galaxy_type'};
		if ( not $type ) {
			warn("No galaxy type for option $name_str given. Skipping");
			return \%result;
		}
		$h{'type'} = $type;

		# extract format
		my $format = $info{'galaxy_format'};
		$h{'format'} = $format if $format;

		# extract value, if given
		my $value = $info{'galaxy_value'};
		$h{'value'} = $value if $value;
		
		$h{'name'} =  $name_str;
		$h{'label'} = $info{'galaxy_label'} || $name_str;		
		$h{'help'} = $description;
		
		# extract options for value, if given
		if ( $info{'galaxy_options'} ) {
			for my $op ( @ {$info{'galaxy_options'}} ) {
				$h{'option'} = [] if not $h{'option'};
				push $h{'option'}, {"value" => $op, "selected" => $value eq $op ? "Yes" : "No" };
			}						
		}	
		
		# add conditionals
		if ( $info{'galaxy_condition'} ) {
			$h{'condition'} = $info{'galaxy_condition'};
		} 

		# add optional flag
		if ( $info{'galaxy_optional'} ) {
			$h{'optional'} = 'true'
		}
		# summarize all fields into return value
		%result = ( $tag => \%h );
	} 
	else {
		return ();
	}
		   	
	return \%result;
}

sub get_conditionals {
	my $in = shift;

	my %conditionals;

	# make hash with params by name
	my %params = map { $_->{'name'}=>$_ } @{$in};

	# loop over parameters and check which ones are 'conditional', meaning their existance
	#  depend on the value of another parameter. If so, the parameter gets nested under
	#  the parameter that is referenced by the condition	
	for my $p ( values %params ) {
		# check if this parameter depends on another
		if ( $p->{'condition'} ) {

			# get the name for the parameter that is referenced
			# NOTE: Only one parameter can be referenced, but
			#  there can be multiple values for the referenced parameter!
			( my $ref_paramname ) = keys ($p->{'condition'});
			( my $ref_p ) = grep { $_->{'name'} eq $ref_paramname } @$in; 

			my @values = map { ref($_) ? @$_ : $_ } values ($p->{'condition'} );
			
			# remove 'condition' tag since it is not a galaxy tag
			delete $p->{'condition'};
			
			$conditionals{$ref_paramname} = {'name'=>"condition_$ref_paramname"} if not $conditionals{$ref_paramname};
			$conditionals{$ref_paramname}->{'param'} = $ref_p if not $conditionals{$ref_paramname}->{'param'};
			$conditionals{$ref_paramname}->{'when'} = [] if not $conditionals{$ref_paramname}->{'when'};

			my %existing = map {$_=>1} map {$_->{'value'}} @{$conditionals{$ref_paramname}->{'when'}};

			# some <when> tags could be already initialized, if another parameter is dependant on 
			#  the reference parameters. Get the values for which <when> tags exist
			for my $v (@values) {
				# make shallow copy of parameter hash
				my $p_add = {%$p};				
				# the <when> tag for this value does not exist yet 
				if ( ! $existing{$v} ) {
					push @{ $conditionals{$ref_paramname}->{'when'} }, {'value' => $v, 'param' => [$p_add]};
				}
				else {
					# iterate over existing <when> tags and add parameter for this value
					for my $w ( @{$conditionals{$ref_paramname}->{'when'}} ) {
						if ( my $val = $w->{'value'} ) {
							if  ( $val eq $v ) {
								# add parameter to existig ones for this value
								push $w->{'param'}, $p_add;
							}
						}
					}
				}
			}	
		}
	}
	return \%conditionals;
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
