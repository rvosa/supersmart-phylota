#!/usr/bin/perl
use strict;
use warnings;
use Template;
use FindBin '$Bin';
use Test::More 'no_plan';
use Bio::Phylo::IO 'parse';

# maybe provided?
my $outfile = shift;

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
	'data'       => $project,
	'ngens'      => 10_000_000,
	'sfreq'      => 10_000,
	'tree_file'  => $tree_file,
	'param_file' => $param_file,
}, \$output ) || die $tt->error();

ok( $output );

if ( $outfile ) {
	open my $fh, '>', $outfile or die $!;
	print $fh $output;
}