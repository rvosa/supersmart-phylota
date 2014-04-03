# This manifest installs the supersmart pipeline, 
# for more information http://www.supersmart-project.org

# update the $PATH environment variable for the Exec tasks.
Exec {
	path => [ 
		"/usr/local/sbin", 
		"/usr/local/bin",
		"/usr/sbin",
		"/usr/bin",
		"/sbin",
		"/bin",
	]
}

# install packages. most likely this is done with yum.
package {
	"mysql-server":             ensure => installed;
	"mysql":                    ensure => installed;
	"ncbi-blast+.x86_64":       ensure => installed;
	"wget":                     ensure => installed;
	"tar":                      ensure => installed;
	"perl-DBI":                 ensure => installed;
	"perl-DBD-MySQL":           ensure => installed;
	"perl-DBIx-Class":          ensure => installed;
	"perl-JSON":                ensure => installed;
	"perl-Moose":               ensure => installed;
	"perl-XML-Twig":            ensure => installed;
	"perl-HTML-Parser":         ensure => installed;
	"perl-Config-Tiny":         ensure => installed;
	"perl-IO-String":           ensure => installed;
	"git":                      ensure => installed;
	"zlib-devel":               ensure => installed;
	"autoconf":                 ensure => installed;
	"automake":                 ensure => installed;
	"libtool":                  ensure => installed;
	"gcc-c++":                  ensure => installed;	
	"curl":                     ensure => installed;
	"gzip":                     ensure => installed;
        "subversion":               ensure => installed;
}

# set up the mysql daemon process
service { 
	"mysqld":
		enable  => true,
		ensure  => running,
		require => Package["mysql-server"],
}

# create links for executables and data directories
file {
	"muscle_link":
		path    => "/usr/local/bin/muscle",
		ensure  => link,
		target  => "/usr/local/src/muscle3.8.31_i86linux64",
		require => Exec["unzip_muscle"];
	"consense_link":
		path    => "/usr/local/bin/consense",
		ensure  => link,
		target  => "/usr/local/src/phylip-3.695/src/consense",
		require => Exec["make_phylip"];
	"examl_link":
		path    => "/usr/local/bin/examl",
		ensure  => link,
		target  => "/usr/local/src/ExaML/examl/examl",
		require => Exec["compile_examl"];
	"parser_link":
		path    => "/usr/local/bin/parser",
		ensure  => link,
		target  => "/usr/local/src/ExaML/parser/parser",
		require => Exec["compile_parser"];
	"treepl_link":
		path    => "/usr/local/bin/treePL",
		ensure  => link,
		target  => "/usr/local/src/treePL/src/treePL",
		require => Exec["compile_treepl"];
	"data_dir":
		path    => "/usr/share/supersmart",
		ensure  => directory;
	"inparanoid_dir":
		path    => "/usr/share/supersmart/inparanoid",
		ensure  => directory;
}

# command line tasks
exec {
  
	# make phylota database
        "dl_phylota_dump":
		command => "wget https://dl.dropboxusercontent.com/u/4180059/phylota.tar.gz",
		cwd     => "/usr/share/supersmart",
		creates => "/usr/share/supersmart/phylota.tar.gz",
		timeout => 0,
		require => [ File[ 'data_dir' ], Package[ 'wget' ] ];
	"unzip_phylota_dump":
 		command => "tar -xzvf phylota.tar.gz",
 		creates => "/usr/share/supersmart/phylota",
 		cwd     => "/usr/share/supersmart/",
                timeout => 0,
 		require => Exec[ 'dl_phylota_dump'];
 	"symlink_phylota_dump":
 	        command => "ln -s /usr/share/supersmart/phylota",
 		creates => "/var/lib/mysql/phylota/seqs.frm",
 		cwd     => "/var/lib/mysql/",
 		require => Exec[ 'unzip_phylota_dump' ];
 	"chown_phylota_db":
 		command => "chown -R -h mysql:mysql phylota/",
 		cwd     => "/var/lib/mysql/",
 		require => Exec[ 'symlink_phylota_dump' ];
 	"grant_phylota_db":
 		command => "mysql -u root -e \"grant all on phylota.* to 'mysql'@'localhost';\"",
 		require => Exec[ 'chown_phylota_db' ];
 		
 	# install mafft
 	"dl_mafft":
 		command => "wget http://mafft.cbrc.jp/alignment/software/mafft-7.130-without-extensions-src.tgz",
 		cwd     => "/usr/share/supersmart",
 		creates => "/usr/share/supersmart/mafft-7.130-without-extensions-src.tgz",
		require => [ Package[ 'wget' ] ]; 		
	"unzip_mafft":
		command => "tar -xzvf mafft-7.130-without-extensions-src.tgz",
 		cwd     => "/usr/share/supersmart",
 		creates => "/usr/share/supersmart/mafft-7.130-without-extensions",		
 		require => Exec[ 'dl_mafft' ];
 	"install_mafft":
 		command => "make install",
 		cwd     => "/usr/share/supersmart/mafft-7.130-without-extensions/core",		
 		require => Exec[ 'unzip_mafft' ],
 		creates => "/usr/local/bin/mafft"; 		

	# make inparanoid blast db
	"dl_inparanoid_seq":
		command => "curl http://inparanoid.sbc.su.se/download/current/sequences/processed/FASTA -o inparanoid.fa",
		cwd     => "/usr/share/supersmart",
		creates => "/usr/share/supersmart/inparanoid.fa",
		timeout => 0,		
		require => [ File[ 'data_dir' ], Package[ 'curl' ] ];	
	"inparanoid_formatdb":
		command => "makeblastdb -in inparanoid.fa -logfile inparanoid-blast.log -dbtype prot -parse_seqids",
		cwd     => "/usr/share/supersmart",
		creates => "/usr/share/supersmart/inparanoid-blast.log",
		timeout => 0,
		require => [ Exec[ 'dl_inparanoid_seq' ], Package[ 'ncbi-blast+.x86_64' ] ];		

	# install muscle multiple sequence alignment
	"download_muscle":
		command => "wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/muscle3.8.31_i86linux64.tar.gz",
		require => Package[ 'wget', 'tar' ];
	"unzip_muscle":
		command => "tar -xzvf /usr/local/src/muscle3.8.31_i86linux64.tar.gz",
		cwd     => "/usr/local/src",		
		creates => "/usr/local/src/muscle3.8.31_i86linux64",
		require => Exec["download_muscle"];

	# install Bio::Phylo
	"clone_bio_phylo":
		command => "git clone https://github.com/rvosa/bio-phylo.git",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/bio-phylo",
		require => Package[ 'git' ];
	"make_bio_phylo_sh":
		command => "echo 'export PERL5LIB=\$PERL5LIB:/usr/local/src/bio-phylo/lib' > biophylo.sh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/biophylo.sh",
		require => Exec[ 'clone_bio_phylo' ];
	"make_bio_phylo_csh":
	  command => "echo 'setenv PERL5LIB \$PERL5LIB:/usr/local/src/bio-phylo/lib' > biophylo.csh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/biophylo.csh",
		require => Exec[ 'clone_bio_phylo' ];
	
	# install bioperl-live
	"clone_bioperl_live":
		command => "git clone -b v1.6.x https://github.com/bioperl/bioperl-live.git",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/bioperl-live",
		require => Package[ 'git' ];		
	"make_bioperl_live_sh":
		command => "echo 'export PERL5LIB=\$PERL5LIB:/usr/local/src/bioperl-live' > bioperl_live.sh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/bioperl_live.sh",
		require => Exec[ 'clone_bioperl_live' ];
	"make_bioperl_live_csh":
		command => "echo 'setenv PERL5LIB \$PERL5LIB:/usr/local/src/bioperl-live' > bioperl_live.csh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/bioperl_live.csh",
		require => Exec[ 'clone_bioperl_live' ];		

	# install bioperl-run
	"clone_bioperl_run":
		command => "git clone https://github.com/bioperl/bioperl-run.git",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/bioperl-run",
		require => Package[ 'git' ];		
	"make_bioperl_run_sh":
		command => "echo 'export PERL5LIB=\$PERL5LIB:/usr/local/src/bioperl-run/lib' > bioperl_run.sh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/bioperl_run.sh",
		require => Exec[ 'clone_bioperl_run' ];
	"make_bioperl_run_csh":
		command => "echo 'setenv PERL5LIB \$PERL5LIB:/usr/local/src/bioperl-run/lib' > bioperl_run.csh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/bioperl_run.csh",
		require => Exec[ 'clone_bioperl_run' ];	
		
        #install perl package Math::Random
        "download_math_random":
                command => "wget http://search.cpan.org/CPAN/authors/id/G/GR/GROMMEL/Math-Random-0.70.tar.gz",
                cwd     => "/usr/local/src",
                creates => "/usr/local/src/Math-Random-0.70.tar.gz",
                require => Package[ 'wget', 'tar' ];
        "unzip_math_random":
                command => "tar -xvzf /usr/local/src/Math-Random-0.70.tar.gz",
                creates => "/usr/local/src/Math-Random-0.70/Makefile.PL",
                cwd     => "/usr/local/src",
                require => Exec["download_math_random"];
        "make_makefile_math_random":
                command => "perl Makefile.PL",
		cwd     => "/usr/local/src/Math-Random-0.70",
		creates => "/usr/local/src/Math-Random-0.70/Makefile",
		require => Exec["unzip_math_random"];
        "make_install_math_random":
		command => "make install LD_LIBRARY_PATH=/usr/lib",
		cwd     => "/usr/local/src/Math-Random-0.70",		
		creates => "/usr/local/lib64/perl5/Math/Random.pm",
		require => Exec["make_makefile_math_random"];	
          
        # install perl package Parallel::MPI::Simple
	"download_parallel_mpi_simple":
		command => "wget http://search.cpan.org/CPAN/authors/id/A/AJ/AJGOUGH/Parallel-MPI-Simple-0.10.tar.gz",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/Parallel-MPI-Simple-0.10.tar.gz",
		require => Package[ 'wget', 'tar' ];		
	"unzip_parallel_mpi_simple":
		command => "tar -xzvf /usr/local/src/Parallel-MPI-Simple-0.10.tar.gz",
		cwd     => "/usr/local/src",		
		creates => "/usr/local/src/Parallel-MPI-Simple-0.10/Makefile.PL",
		require => Exec["download_parallel_mpi_simple"];
	"make_makefile_parallel_mpi_simple":
		command => "perl Makefile.PL CCFLAGS='-I/usr/include' LDFLAGS='-L/usr/lib/ -lmpi' CC=mpicc LD=mpicc",
		cwd     => "/usr/local/src/Parallel-MPI-Simple-0.10/",
		creates => "/usr/local/src/Parallel-MPI-Simple-0.10/Makefile",
		require => Exec["unzip_parallel_mpi_simple","install_openmpi"];
	"make_install_parallel_mpi_simple":
		command => "make install LD_LIBRARY_PATH=/usr/lib",
		cwd     => "/usr/local/src/Parallel-MPI-Simple-0.10/",		
		creates => "/usr/local/src/Parallel-MPI-Simple-0.10/Simple.so",
		require => Exec["make_makefile_parallel_mpi_simple"];	
		
	# install openmpi
	"download_openmpi":
		command => "wget http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.5.tar.gz",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/openmpi-1.6.5.tar.gz",
		require => Package[ 'wget', 'tar' ];
	"unzip_openmpi":
		command => "tar -xvzf openmpi-1.6.5.tar.gz",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/openmpi-1.6.5",
		require => Exec['download_openmpi'];
	"install_openmpi":
		command => "/usr/local/src/openmpi-1.6.5/configure --prefix=/usr --disable-dlopen && make install",
		creates => "/usr/bin/mpicc",
		cwd     => "/usr/local/src/openmpi-1.6.5",
		timeout => 0,
		require => Exec['unzip_openmpi'];
	
	# install phylip
	"download_phylip":
		command => "wget http://evolution.gs.washington.edu/phylip/download/phylip-3.695.tar.gz",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/phylip-3.695.tar.gz",
		require => Package[ 'wget', 'tar' ];
	"unzip_phylip":
		command => "tar -xzvf /usr/local/src/phylip-3.695.tar.gz",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/phylip-3.695/src/Makefile.unx",
		require => Exec["download_phylip"];
	"make_phylip":
		command => "make -f Makefile.unx consense",
		cwd     => "/usr/local/src/phylip-3.695/src/",
		creates => "/usr/local/src/phylip-3.695/src/consense",
		require => Exec["unzip_phylip"];
	
	# install examl & parser
	"clone_examl":
		command => "git clone https://github.com/stamatak/ExaML.git",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/ExaML/",
		require => Package[ 'git', 'zlib-devel' ];
	"compile_examl":
		command => "make -f Makefile.SSE3.gcc LD_LIBRARY_PATH=/usr/lib",
		cwd     => "/usr/local/src/ExaML/examl",
		creates => "/usr/local/src/ExaML/examl/examl",
		require => Exec["clone_examl","install_openmpi"];
	"compile_parser":
		command => "make -f Makefile.SSE3.gcc",
		cwd     => "/usr/local/src/ExaML/parser",
		creates => "/usr/local/src/ExaML/parser/parser",
		require => Exec["clone_examl","install_openmpi"];
		
	# install supersmart
	"clone_supersmart":
		command => "git clone https://github.com/naturalis/supersmart.git",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/supersmart",
		require => Package[ 'git' ];
	"make_supersmart_sh":
		command => "echo 'export LD_LIBRARY_PATH=/usr/lib' > supersmart.sh && echo 'export SUPERSMART_HOME=/usr/local/src/supersmart' >> supersmart.sh && echo 'export PERL5LIB=\$PERL5LIB:\$SUPERSMART_HOME/lib' >> supersmart.sh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/supersmart.sh";
	"make_supersmart_csh":
		command => "echo 'setenv LD_LIBRARY_PATH /usr/lib' > supersmart.csh && echo 'setenv SUPERSMART_HOME /usr/local/src/supersmart' >> supersmart.csh && echo 'setenv PERL5LIB \$PERL5LIB:\$SUPERSMART_HOME/lib' >> supersmart.csh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/supersmart.csh";		
        "chown_supersmart_home":
                command => "chown -R vagrant .",
                cwd     => "/usr/local/src/supersmart",
                require => Exec["clone_supersmart"];

        # install BEAST
        "download_beast":
                command => "wget https://beast-mcmc.googlecode.com/files/BEASTv1.8.0.tgz",
                cwd     => "/usr/local/src",
                creates => "/usr/local/src/BEASTv1.8.0.tgz",
                require => Package[ 'wget' ];
        "unzip_beast":
                command => "tar -xvzf BEASTv1.8.0.tgz",
                cwd     => "/usr/local/src/",
                creates => "/usr/local/src/BEASTv1.8.0/bin/beast",
                require => Exec[ 'download_beast' ];

        # install beagle-lib
        "checkout_beagle_lib":
                command => "svn checkout http://beagle-lib.googlecode.com/svn/trunk/ beagle-lib",
                cwd     => "/usr/local/src",
                creates => "/usr/local/src/beagle-lib/autogen.sh",
                require => Package[ 'subversion' ];
        "generate_beagle_config":
                command => "sh autogen.sh",
                cwd     => "/usr/local/src/beagle-lib/",
                creates => "/usr/local/src/beagle-lib/configure",
                require => Exec[ 'checkout_beagle_lib' ];
        "build_beagle_lib":
                command => "sh configure && make install",
                cwd     => "/usr/local/src/beagle-lib/",
                creates => "/usr/local/lib/libhmsbeagle.so",
                require => Exec[ 'generate_beagle_config' ];
          
	# install treePL
	"clone_treepl":
		command => "git clone https://github.com/blackrim/treePL.git",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/treePL",
		require => Package[ 'git', 'autoconf', 'automake', 'libtool', 'gcc-c++', 'tar' ];
	"unzip_adolc":
		command => "tar -xvzf adol-c_git_saved.tar.gz",
		cwd     => "/usr/local/src/treePL/deps",
		creates => "/usr/local/src/treePL/deps/adol-c",
		require => Exec["clone_treepl","install_openmpi"];
	"autoreconf_adolc":
		command => "autoreconf -fi",
		cwd     => "/usr/local/src/treePL/deps/adol-c",
		creates => "/usr/local/src/treePL/deps/adol-c/configure",
		require => Exec["unzip_adolc"];
	"install_adolc":
		command => "/usr/local/src/treePL/deps/adol-c/configure --prefix=/usr --with-openmp-flag=-fopenmp && make install",
		cwd     => "/usr/local/src/treePL/deps/adol-c",
		creates => "/usr/lib64/libadolc.so",
		require => Exec["autoreconf_adolc"];
	"unzip_nlopt":
		command => "tar -xzvf nlopt-2.2.4.tar.gz",
		cwd     => "/usr/local/src/treePL/deps",
		creates => "/usr/local/src/treePL/deps/nlopt-2.2.4",
		require => Exec["clone_treepl","install_openmpi"];
	"install_nlopt":
		command => "/usr/local/src/treePL/deps/nlopt-2.2.4/configure && make install",
		cwd		=> "/usr/local/src/treePL/deps/nlopt-2.2.4",
		creates => "/usr/local/include/nlopt.h",
		require => Exec["unzip_nlopt"];
	"compile_treepl":
		command => "/usr/local/src/treePL/src/configure && make",
		cwd     => "/usr/local/src/treePL/src",
		creates => "/usr/local/src/treePL/src/treePL",
		require => Exec["install_adolc","install_nlopt"];		
}
