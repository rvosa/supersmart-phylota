PRE-REQUISITES

The pipeline depends on a number of additional CPAN modules. These are listed in
http://search.cpan.org/dist/Bundle-Bio-Phylo-PhyLoTA/, so that they can all (in
principle) be installed using:

$ perl -MCPAN -e 'install Bundle::Bio::Phylo::PhyLoTA'

That said, some of these installs (e.g. bioperl-live, bioperl-run, DBD::mysql)
may be more complicated than that. Consult their respective documentation or
seek help from your sysadmin.

In addition to these CPAN modules, there's a number of compiled 3rd party
executables that either need to be in your system's PATH or specified in the
config file (see below):

- muscle
- examl
- parser
- phytime
- formatdb
- blastall
- consense (from PHYLIP)

CONFIGURATION

The pipeline API has a single point of configuration: the conf/phylota.ini file.
Locations of 3rd party executables, run-time parameters, input and output file
names etc. are defined there (and only there).

PARALLELIZATION

At present, the pipeline is designed to be parallelized on an MPI architecture,
specifically OpenMPI.
