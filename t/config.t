#!/usr/bin/perl
use Test::More 'no_plan';
use Net::Ping;
use Config;
use File::Spec;

BEGIN { use_ok('Bio::Phylo::PhyLoTA::Config'); }
my $config = new_ok('Bio::Phylo::PhyLoTA::Config');

# verify these files exist
my @files = qw(
	GB_RELNUM_FILE
);
for my $file ( @files ) {
	ok( $config->$file, "file $file is defined" );
	ok( -f $config->$file, "file $file exists" );
}

# verify these binaries can be found
my @bins = qw(
	EXAML_BIN
	EXABAYES_PARSER_BIN
        EXAML_PARSER_BIN
	BLASTP_BIN
	MAKEBLASTDB_BIN
	MUSCLE_BIN
	PHYLIP_CONSENSE_BIN
	EXABAYES_CONSENSE_BIN
        TREEPL_BIN
);
BIN: for my $var ( @bins ) {
	my $bin = $config->$var;
	if ( -x $bin ) {
		ok( -x $bin, "$bin is a path" );
	}
	else {
		for my $dir ( split /$Config{path_sep}/, $ENV{'PATH'} ) {
			my $file = File::Spec->catfile( $dir, $bin . $Config{'_exe'} );
			if ( -x $file ) {
				ok( -x $file, "$file found on PATH" );
				next BIN;
			}
			else {
				#warn "$file doesn't exist";
			}
		}
		ok( 0, "$bin not found" ); # will evaluate to failure
	}
}


# verify these servers are reachable
my @servers = qw(HOST SERVER);
my $p = Net::Ping->new;
for my $server ( @servers ) {
	SKIP : {
		skip "no SERVER available yet", 2 if $server =~ /SERVER/;
		ok( $config->$server, "server $server is defined" );
		ok( $p->ping($config->$server), "server $server is reachable" );
	};
}
