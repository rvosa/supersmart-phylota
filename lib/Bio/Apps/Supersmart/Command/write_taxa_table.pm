package Bio::Apps::Supersmart::Command::write_taxa_table;

# ABSTRACT: Writes taxa table for given list of taxon names or root taxa.

use Bio::Apps::Supersmart -command;


sub opt_spec {
	return (
		["infile|i", "file with list of taxon names", ],
		["expand_rank|e", "taxonomic rank to which root taxa are expanded", { default => 0}],	
	);
	
}


sub execute {
	print "executing write_taxa_table\n";		
}



1;
