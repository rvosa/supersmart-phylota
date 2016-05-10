#!/usr/bin/perl
use strict;
use warnings;
use Template;
use FindBin '$Bin';
use Test::More 'no_plan';
use Bio::Phylo::IO 'parse';

# maybe provided?
my $beastfile_name = shift;

# test file locations
my $template = $ENV{'SUPERSMART_HOME'} . '/data/BEAST/starbeast.xml';
my $data = $Bin . '/testdata/clade1.xml';
ok( -e $template, "found $template" );
ok( -e $data, "found $data" );

# create outfile locations
my $stem = $data;
$stem =~ s/\..+$//;
my $tree_file  = $stem . '.nex';
my $param_file = $stem . '.log';

# read prepared NeXML
my $project = parse(
	'-format'     => 'nexml',
	'-file'       => $data,
	'-as_project' => 1,
);
isa_ok( $project, 'Bio::Phylo::Project' );

# instantiate Template::Toolkit
my $tt = Template->new({ 'ABSOLUTE' => 1 });
isa_ok( $tt, 'Template' );

my $output;
$tt->process( $template, {
	'data'         => $project,
	'chain_length' => 10_000_000,
	'sample_freq'  => 10_000,
	'outfile_name' => $tree_file,
	'logfile_name' => $param_file,
}, \$output ) || die $tt->error();

ok( $output );

if ( $beastfile_name ) {
	open my $fh, '>', $beastfile_name or die $!;
	print $fh $output;
}