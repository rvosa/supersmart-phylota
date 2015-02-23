#!/usr/bin/perl
use strict;
use warnings;
use Test::More;
use Bio::Phylo::PhyLoTA::Service;

eval "use Test::Pod::Coverage 1.00";
plan skip_all => "Test::Pod::Coverage 1.00 required for testing POD coverage" if $@;

# this will prevent pod-coverage from complaining about missing documentation for subroutines that are 
#  overridden in child classes
my $inherited = { trustme => [qr/^(options|run|validate)$/] };

all_pod_coverage_ok($inherited);