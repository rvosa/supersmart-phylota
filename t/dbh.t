#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';

# this is a unit test
BEGIN { use_ok('Bio::Phylo::PhyLoTA::DBH'); }
my $dbh = new_ok('Bio::Phylo::PhyLoTA::DBH');

ok( $dbh->ping, "connected to database" );

for my $field ( qw(host dsn database rdbms user) ) {
	ok( $dbh->$field, "$field is defined" );
}


