#!/usr/bin/perl
use strict;
use warnings;
use lib 'lib';
use ExtUtils::MakeMaker;

WriteMakefile(
	'NAME' => 'Bio::Phylo::PhyLoTA',
	'PREREQ_PM' => {
		'DBIx::Class'           => 0,
		'DBI'                   => 0,
		'DBD::mysql'            => 0,
		'Template'              => 0,
		'Config::Tiny'          => 0,
		'Moose'                 => 0,
		'XML::Twig'             => 0,
		'HTML::Parser'          => 0,
		'JSON'                  => 0,
		'Math::Random'          => 0,
		'App::Cmd'              => 0,
		'String::RewritePrefix' => 0,
		'IO::String'            => 0,
		'App::Cmd'              => 0,
		# 'BioPerl'               => 0,
		# 'BioPerl-Run'           => 0,
		# 'Parallel::MPI::Simple' => 0, 
		# 'Bio::Phylo'            => 0, 
	}
);

