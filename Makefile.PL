#!/usr/bin/perl
use strict;
use warnings;
use lib 'lib';
use ExtUtils::MakeMaker;

WriteMakefile(
	'NAME' => 'Bio::Phylo::PhyLoTA',
	'PREREQ_PM' => {
		'DBIx::Class'           => 0,
		'Getopt::Long'          => 0,
		'Pod::Usage'            => 0,
		'Bio::Phylo'            => 0,
		'Template'              => 0,
		'Config::Tiny'          => 0,
		'DBI'                   => 0,
		'DBD::mysql'            => 0,
		'Moose'                 => 0,
		'Bio::Species'          => 0,
		'XML::Twig'             => 0,
		'HTML::Parser'          => 0,
		'JSON'                  => 0,
		'Parallel::MPI::Simple' => 0,
	}
);

