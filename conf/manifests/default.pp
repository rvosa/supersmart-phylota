# This manifest installs the supersmart pipeline, 
# for more information http://www.supersmart-project.org

# set user and default paths for storing data, tools and source code
$username = "vagrant"
$supersmart_dir	= "/home/${username}/SUPERSMART"
$tools_dir		= "${supersmart_dir}/tools"
$tools_bin_dir	= "${tools_dir}/bin"
$src_dir		= "${supersmart_dir}/src"
$data_dir		= "${supersmart_dir}/data"

# update the $PATH environment variable
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

# install packages.
package {
	"mysql-server":                ensure => installed;
	"mysql-client":                ensure => installed;
	"ncbi-blast+":       	       ensure => installed;
	"wget":                        ensure => installed;
	"tar":                         ensure => installed;
	"libdbi-perl":                 ensure => installed;
	"libdbd-mysql-perl":           ensure => installed;
	"libdbix-class-perl":          ensure => installed;
	"libjson-perl":                ensure => installed;
	"libmoose-perl":               ensure => installed;
	"libxml-twig-perl":            ensure => installed;
	"libhtml-html5-parser-perl":   ensure => installed;
	"libconfig-tiny-perl":         ensure => installed;
	"libio-string-perl":           ensure => installed;
	"git":                         ensure => installed;
	"libarchive-dev":              ensure => installed;
  	"zlib1g-dev":                  ensure => installed;
 	"autoconf":                    ensure => installed;
	"automake":                    ensure => installed;
	"libtool":                     ensure => installed;
	"build-essential":             ensure => installed;	
	"curl":                        ensure => installed;
	"gzip":                        ensure => installed;
  	"openjdk-6-jdk":               ensure => installed;
 	"subversion":                  ensure => installed;
}

# set up the mysql daemon process
service { 
	"mysql":
		enable  => true,
		ensure  => running,
		require => Package["mysql-server"],
}

# create links for executables and data directories
file {
	$supersmart_dir:
		ensure  => directory,
		group   => $username,
		owner   => $username;
	$data_dir:
		ensure  => directory,
    	group   => $username,
    	owner   => $username;
	$src_dir:
		ensure  => directory,
    	group   => $username,
    	owner   => $username,
    	recurse => true;
	$tools_dir:
		ensure  => directory,
	  	group   => $username,
    	owner   => $username;
  	$tools_bin_dir:
		ensure  => directory,
	  	group   => $username,
    	owner   => $username;

  "muscle_link":
		path    => "${tools_bin_dir}/muscle",
		ensure  => link,
		target  => "${tools_bin_dir}/muscle3.8.31_i86linux64",
		require => Exec["unzip_muscle"];
	"phylip-consense_link":
		path    => "/usr/local/bin/consense-phylip",
		ensure  => link,
		target  => "${tools_dir}/phylip-3.696/src/consense",
		require => Exec["make_phylip"];
	"examl_link":
		path    => "${tools_bin_dir}/examl",
		ensure  => link,
		target  => "${tools_dir}/ExaML/examl/examl",
		require => Exec["compile_examl"];
  	"parse_examl_link":
    	path    => "${tools_bin_dir}/parse-examl",
    	ensure  => link,
    	target  => "${tools_dir}/ExaML/parser/parse-examl",
    	require => Exec["compile_examl"];     
	"exabayes_link":
		path    => "${tools_bin_dir}/exabayes",
		ensure  => link,
		target  => "${tools_dir}/exabayes-1.2.1/bin/exabayes",
		require => Exec["compile_exabayes"];
	"parse_exabayes_link":
    	path  => "${tools_bin_dir}/parse-exabayes",
    	ensure  => link,
    	target  => "${tools_dir}/exabayes-1.2.1/bin/parser",
    	require => Exec["compile_exabayes"];     	 
   "exabayes_consense_link":
    	path  => "${tools_bin_dir}/consense-exabayes",
    	ensure  => link,
    	target  => "${tools_dir}/exabayes-1.2.1/bin/consense",
    	require => Exec["compile_exabayes"];       
   "raxml_link":
    	path    => "${tools_bin_dir}/raxml",
    	ensure  => link,
    	target  => "${tools_dir}/standard-RAxML/raxmlHPC-AVX",
    	require => Exec["compile_raxml"];		 
   "treepl_link":
	  	path    => "${tools_bin_dir}/treePL",
		ensure  => link,
		target  => "${tools_dir}/treePL/src/treePL",
		require => Exec["compile_treepl"];
}

# command line tasks
exec {
  
   	# add bin directory for all required tools and smrt executables to PATH
  	"make_bindir_sh":
    	command => "echo 'export PATH=\$PATH:${tools_bin_dir}:${src_dir}/supersmart/app' > supersmart-tools-bin.sh",
    	cwd     => "/etc/profile.d",
    	creates => "/etc/profile.d/supersmart-tools-bin.sh",
    	require => File[ $tools_bin_dir ];
  	"make_bindir_csh":
    	command => "echo 'setenv PATH \$PATH:${tools_bin_dir}:${src_dir}/supersmart/app' > supersmart-tools-bin.csh",
    	cwd     => "/etc/profile.d",
    	creates => "/etc/profile.d/supersmart-tools-bin.csh",
    	require => File[ $tools_bin_dir ];
  
  	# add default EDITOR environment variable
  	"make_editor_sh":
    	command => "echo 'export EDITOR=/etc/alternatives/editor' > supersmart-editor.sh",
    	cwd     => "/etc/profile.d",
    	creates => "/etc/profile.d/supersmart-editor.sh";
  	"make_editor_csh":
    	command => "echo 'setenv EDITOR=/etc/alternatives/editor' > supersmart-editor.csh",
    	cwd     => "/etc/profile.d",
    	creates => "/etc/profile.d/supersmart-editor.csh";
  
	# make phylota database
	"dl_phylota_dump":
		command => "wget https://dl.dropboxusercontent.com/u/4180059/phylota.tar.gz",
		cwd     => $data_dir,
		creates => "${data_dir}/phylota.tar.gz",
		timeout => 0,
		require => [ File[ $data_dir ], Package[ 'wget' ] ];
	"unzip_phylota_dump":
 		command => "tar -xzvf phylota.tar.gz",
 		creates => "${data_dir}/phylota",
 		cwd     => $data_dir,
    timeout => 0,
 		require => Exec[ 'dl_phylota_dump'];
 	"mv_phylota_dump":
 	    command => "mv ${data_dir}/phylota/ /var/lib/mysql",
 		creates => "/var/lib/mysql/phylota/seqs.frm",
 		cwd     => "/var/lib/mysql/",
 		require => Exec[ 'unzip_phylota_dump' ];
 	"chown_phylota_db":
 		command => "chown -R -h mysql:mysql phylota/",
 		cwd     => "/var/lib/mysql/",
 		require => Exec[ 'mv_phylota_dump' ];
 	"grant_phylota_db":
 		command => "mysql -u root -e \"grant all on phylota.* to 'mysql'@'localhost';\"",
 		require => Exec[ 'chown_phylota_db' ];
 		
 	# install mafft
 	"dl_mafft":
 		command => "wget http://mafft.cbrc.jp/alignment/software/mafft-7.130-without-extensions-src.tgz",
 		cwd     => $tools_dir,
 		creates => "${tools_dir}/mafft-7.130-without-extensions-src.tgz",
		require => [ File [ $tools_dir ] ]; 		
	"unzip_mafft":
		command => "tar -xzvf mafft-7.130-without-extensions-src.tgz",
 		cwd     => $tools_dir,
 		creates => "${tools_dir}/mafft-7.130-without-extensions",		
 		require => Exec[ 'dl_mafft' ];
 	"make_nonroot_makefile_mafft":
 		command => "sed -i -e 's|PREFIX = /usr/local|PREFIX = ${tools_dir}|g' Makefile",
 		cwd	=> "${tools_dir}/mafft-7.130-without-extensions/core",
 		require => Exec[ 'unzip_mafft' ]; 		 	
 	"install_mafft":
 		command => "make && make install",
 		cwd     => "${tools_dir}/mafft-7.130-without-extensions/core",		
 		creates => "$tools_dir/bin/mafft", 		
 		require => Exec[ 'make_nonroot_makefile_mafft' ];

	# install muscle multiple sequence alignment
	"download_muscle":
		command => "wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz",
		cwd     => $tools_dir,
		creates => "${tools_dir}/muscle3.8.31_i86linux64.tar.gz",
		require => Package[ 'wget', 'tar' ];
	"unzip_muscle":
		command => "tar -xzvf muscle3.8.31_i86linux64.tar.gz -C ${tools_bin_dir}",
		cwd     => $tools_dir,		
		creates => "${tools_bin_dir}/muscle3.8.31_i86linux64",
		require => Exec["download_muscle"];

 	# install Bio::Phylo
	"clone_bio_phylo":
		command => "git clone https://github.com/rvosa/bio-phylo.git",
		cwd     => $tools_dir,
		creates => "${tools_dir}/bio-phylo",
		timeout => 0,
		require => Package[ 'git' ];
	"make_bio_phylo_sh":
		command => "echo 'export PERL5LIB=\$PERL5LIB:${tools_dir}/bio-phylo/lib' > biophylo.sh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/biophylo.sh",
		require => Exec[ 'clone_bio_phylo' ];
	"make_bio_phylo_csh":
	  command => "echo 'setenv PERL5LIB \$PERL5LIB:${tools_dir}/bio-phylo/lib' > biophylo.csh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/biophylo.csh",
		require => Exec[ 'clone_bio_phylo' ];
	
	# install bioperl-live
	"clone_bioperl_live":
		command => "git clone -b v1.6.x https://github.com/bioperl/bioperl-live.git",
		cwd     => $tools_dir,
		creates => "${tools_dir}/bioperl-live",
		timeout => 0,
		require => Package[ 'git' ];		
	"make_bioperl_live_sh":
		command => "echo 'export PERL5LIB=\$PERL5LIB:${tools_dir}/bioperl-live' > bioperl_live.sh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/bioperl_live.sh",
		require => Exec[ 'clone_bioperl_live' ];
	"make_bioperl_live_csh":
		command => "echo 'setenv PERL5LIB \$PERL5LIB:${tools_dir}/bioperl-live' > bioperl_live.csh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/bioperl_live.csh",
		require => Exec[ 'clone_bioperl_live' ];		

	# install bioperl-run
	"clone_bioperl_run":
		# NOTE: using forked version from hettling's repository in which RaXML wrapper works
		#       This can be changed to the official bioperl repo if they accept the changes.
		command => "git clone https://github.com/hettling/bioperl-run.git",
		cwd     => $tools_dir,
		creates => "${tools_dir}/bioperl-run",
		require => Package[ 'git' ];		
	"make_bioperl_run_sh":
		command => "echo 'export PERL5LIB=\$PERL5LIB:${tools_dir}/bioperl-run/lib' > bioperl_run.sh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/bioperl_run.sh",
		require => Exec[ 'clone_bioperl_run' ];
	"make_bioperl_run_csh":
		command => "echo 'setenv PERL5LIB \$PERL5LIB:${tools_dir}/bioperl-run/lib' > bioperl_run.csh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/bioperl_run.csh",
		require => Exec[ 'clone_bioperl_run' ];	
		
  	#install perl package App:Cmd     
   	"clone_app_cmd":
   		command => "git clone https://github.com/rjbs/App-Cmd.git",
   		cwd     => $tools_dir,
   		creates => "${tools_dir}/App-Cmd",
   		require => Package[ 'git' ];  
  	"make_app_cmd_run_sh":
    	command => "echo 'export PERL5LIB=\$PERL5LIB:${tools_dir}/App-Cmd/lib' > app_cmd.sh",
    	cwd     => "/etc/profile.d",
    	creates => "/etc/profile.d/app_cmd.sh",
    	require => Exec[ 'clone_app_cmd' ];
  	"make_app_cmd_run_csh":
    	command => "echo 'setenv PERL5LIB \$PERL5LIB:${tools_dir}/App-Cmd/lib' > app_cmd.csh",
    	cwd     => "/etc/profile.d",
    	creates => "/etc/profile.d/app_cmd.csh",
    	require => Exec[ 'clone_app_cmd' ]; 
  	  
  	#install perl package String::RewritePrefix, needed by App::Cmd     
  	"download_string_rewrite_prefix":
    	command => "wget http://search.cpan.org/CPAN/authors/id/R/RJ/RJBS/String-RewritePrefix-0.006.tar.gz",
    	cwd     => $tools_dir,
    	creates => "${tools_dir}/String-RewritePrefix-0.006.tar.gz",
    	require => Package[ 'wget', 'tar' ];
   	"unzip_string_rewrite_prefix":
    	command => "tar -xvzf ${tools_dir}/String-RewritePrefix-0.006.tar.gz",
    	creates => "${tools_dir}/String-RewritePrefix-0.006/Makefile.PL",
    	cwd     => $tools_dir,
    	require => Exec["download_string_rewrite_prefix"];
  	"make_string_rewrite_prefix_run_sh":
    	command => "echo 'export PERL5LIB=\$PERL5LIB:${tools_dir}/String-RewritePrefix-0.006/lib' > string_rewrite_prefix.sh",
    	cwd     => "/etc/profile.d",
    	creates => "/etc/profile.d/string_rewrite_prefix.sh",
    	require => Exec[ 'unzip_string_rewrite_prefix' ];
  	"make_string_rewrite_prefix_run_csh":
    	command => "echo 'setenv PERL5LIB \$PERL5LIB:${tools_dir}/String-RewritePrefix-0.006/lib' > string_rewrite_prefix.csh",
    	cwd     => "/etc/profile.d",
    	creates => "/etc/profile.d/string_rewrite_prefix.csh",
    	require => Exec[ 'unzip_string_rewrite_prefix' ]; 	
  	  
   	#install perl package Math::Random
  	"download_math_random":
    	command => "wget http://search.cpan.org/CPAN/authors/id/G/GR/GROMMEL/Math-Random-0.70.tar.gz",
    	cwd     => $tools_dir,
    	creates => "${tools_dir}/Math-Random-0.70.tar.gz",
    	require => Package[ 'wget', 'tar' ];
   	"unzip_math_random":
    	command => "tar -xvzf ${tools_dir}/Math-Random-0.70.tar.gz",
    	creates => "${tools_dir}/Math-Random-0.70/Makefile.PL",
    	cwd     => $tools_dir,
    	require => Exec["download_math_random"];
   "make_makefile_math_random":
    	command => "perl Makefile.PL",
		cwd     => "${tools_dir}/Math-Random-0.70",
		creates => "${tools_dir}/Math-Random-0.70/Makefile",
		require => Exec["unzip_math_random"];
   "make_install_math_random":
		command => "make install LD_LIBRARY_PATH=/usr/lib",
		cwd     => "${tools_dir}/Math-Random-0.70",		
		creates => "/usr/local/lib64/perl5/Math/Random.pm",
		require => Exec["make_makefile_math_random"];	
                                   	      
	# install openmpi
	"download_openmpi":
		command => "wget http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.5.tar.gz",
		cwd     => $tools_dir,
		creates => "${tools_dir}/openmpi-1.6.5.tar.gz",
		require => Package[ 'wget', 'tar' ];
	"unzip_openmpi":
		command => "tar -xvzf openmpi-1.6.5.tar.gz",
		cwd     => $tools_dir,
		creates => "${tools_dir}/openmpi-1.6.5",
		require => Exec['download_openmpi'];
	"install_openmpi":
		command => "${tools_dir}/openmpi-1.6.5/configure --prefix=/usr --disable-dlopen && make install",
		creates => "/usr/bin/mpicc",
		cwd     => "${tools_dir}/openmpi-1.6.5",
		timeout => 0,
		require => Exec['unzip_openmpi'];
	
	# install phylip
	"download_phylip":
		command => "wget http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz",
		cwd     => $tools_dir,
		creates => "${tools_dir}/phylip-3.696.tar.gz",
		require => Package[ 'wget', 'tar' ];
	"unzip_phylip":
		command => "tar -xzvf ${tools_dir}/phylip-3.696.tar.gz",
		cwd     => $tools_dir,
		creates => "${tools_dir}/phylip-3.696/src/Makefile.unx",
		require => Exec["download_phylip"];
	"make_phylip":
		command => "make -f Makefile.unx consense",
		cwd     => "${tools_dir}/phylip-3.696/src/",
		creates => "${tools_dir}/phylip-3.696/src/consense",
		require => Exec["unzip_phylip"];
	
	# install examl & parser
	"clone_examl":
		command => "git clone https://github.com/stamatak/ExaML.git",
		cwd     => $tools_dir,
		creates => "${tools_dir}/ExaML/",
	  require => Package[ 'git', 'zlib1g-dev' ];
	"compile_examl":
		command => "make -f Makefile.SSE3.gcc LD_LIBRARY_PATH=/usr/lib",
		cwd     => "${tools_dir}/ExaML/examl",
		creates => "${tools_dir}/ExaML/examl/examl",
		require => Exec["clone_examl","install_openmpi"];
	"compile_parser":
		command => "make -f Makefile.SSE3.gcc",
		cwd     => "${tools_dir}/ExaML/parser",
		creates => "${tools_dir}/ExaML/parser/parser",
		require => Exec["clone_examl","install_openmpi"];

	# install raxml
  	"clone_raxml":
    	command => "git clone https://github.com/stamatak/standard-RAxML.git",
    	cwd     => $tools_dir,
    	creates => "${tools_dir}/standard-RAxML/",
        require => Package[ 'git' ];
  	"compile_raxml":
    	command => "make -f Makefile.AVX.gcc LD_LIBRARY_PATH=/usr/lib",
    	cwd     => "${tools_dir}/standard-RAxML/",
    	creates => "${tools_dir}/standard-RAxML/raxmlHPC-AVX",
    	require => Exec["clone_raxml"];
	 
		  
    # install exabayes
   	"download_exabayes":
    	command => "wget http://sco.h-its.org/exelixis/material/exabayes/1.2.1/exabayes-1.2.1-linux-openmpi-avx.tar.gz",
     	cwd     => $tools_dir,
     	creates => "${tools_dir}/exabayes-1.2.1-linux-openmpi-avx.tar.gz",
     	require => Exec[ "install_openmpi" ];
    "unzip_exabayes":
     	command => "tar -xvzf exabayes-1.2.1-linux-openmpi-avx.tar.gz",
     	cwd     => $tools_dir,
     	creates => "${tools_dir}/exabayes-1.2.1/",
     	require => Exec[ "download_exabayes" ];
    "compile_exabayes":
      	command => "sh build.sh CC=mpicc CXX=mpicxx",
      	cwd     => "${tools_dir}/exabayes-1.2.1/",
      	creates => "${tools_dir}/exabayes-1.2.1/bin/exabayes",
      	timeout => 0,
      	require => Exec[ "unzip_exabayes" ];      
	      
	# install supersmart
	"clone_supersmart":
		command => "git clone https://github.com/naturalis/supersmart.git",
		cwd     => $src_dir,
		creates => "${src_dir}/supersmart",
		user	=> $username,
		require => Package[ 'git' ];
	"make_supersmart_sh":
		command => "echo 'export LD_LIBRARY_PATH=/usr/lib:/usr/lib64:/usr/local/lib' > supersmart.sh && echo 'export SUPERSMART_HOME=${src_dir}/supersmart' >> supersmart.sh && echo 'export PERL5LIB=\$PERL5LIB:\$SUPERSMART_HOME/lib' >> supersmart.sh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/supersmart.sh";
	"make_supersmart_csh":
		command => "echo 'setenv LD_LIBRARY_PATH /usr/lib:/usr/lib64:/usr/local/lib' > supersmart.csh && echo 'setenv SUPERSMART_HOME ${src_dir}/supersmart' >> supersmart.csh && echo 'setenv PERL5LIB \$PERL5LIB:\$SUPERSMART_HOME/lib' >> supersmart.csh",
		cwd     => "/etc/profile.d",
		creates => "/etc/profile.d/supersmart.csh";		
    		    
  	# install BEAST
  	"download_beast":
    	command => "wget https://beast-mcmc.googlecode.com/files/BEASTv1.8.0.tgz",
    	cwd     => $tools_dir,
    	creates => "${tools_dir}/BEASTv1.8.0.tgz",
    	require => Package[ 'wget' ];
  	"unzip_beast":
    	command => "tar -xvzf BEASTv1.8.0.tgz",
    	cwd     => "$tools_dir",
    	creates => "${tools_dir}/BEASTv1.8.0/bin/beast",
    	require => Exec[ 'download_beast' ];

  	"symlink_beast":
    	command => "ln -s ${tools_dir}/BEASTv1.8.0/bin/* .",
    	cwd     => "/usr/local/bin/",
    	creates => "/usr/local/bin/beast",
    	require => Exec[ 'unzip_beast' ];
        
  	# install beagle-lib
  	"checkout_beagle_lib":
    	command => "svn checkout http://beagle-lib.googlecode.com/svn/trunk/ beagle-lib",
    	cwd     => $tools_dir,
    	creates => "${tools_dir}/beagle-lib/autogen.sh",
    	require => Package[ 'subversion', "openjdk-6-jdk" ];
  	"generate_beagle_config":
   		command => "sh autogen.sh",
    	cwd     => "${tools_dir}/beagle-lib/",
    	creates => "${tools_dir}/beagle-lib/configure",
    	require => Exec[ 'checkout_beagle_lib' ];
  	"build_beagle_lib":
    	command => "sh configure && make install",
    	cwd     => "${tools_dir}/beagle-lib/",
    	creates => "/usr/local/lib/libhmsbeagle.so",
    	require => Exec[ 'generate_beagle_config' ];
           
	# install treePL
	"clone_treepl":
		command => "git clone https://github.com/blackrim/treePL.git",
		cwd     => $tools_dir,
		creates => "${tools_dir}/treePL",
		require => Package[ 'git', 'autoconf', 'automake', 'libtool', 'build-essential', 'tar' ];
	"unzip_adolc":
		command => "tar -xvzf adol-c_git_saved.tar.gz",
		cwd     => "${tools_dir}/treePL/deps",
		creates => "${tools_dir}/treePL/deps/adol-c",
		require => Exec["clone_treepl","install_openmpi"];
	"autoreconf_adolc":
		command => "autoreconf -fi",
		cwd     => "${tools_dir}/treePL/deps/adol-c",
		creates => "${tools_dir}/treePL/deps/adol-c/configure",
		require => Exec["unzip_adolc"];
	"install_adolc":
		command => "${tools_dir}/treePL/deps/adol-c/configure --prefix=/usr --with-openmp-flag=-fopenmp && make install",
		cwd     => "${tools_dir}/treePL/deps/adol-c",
		creates => "/usr/lib64/libadolc.so",
		require => Exec["autoreconf_adolc"];
	"unzip_nlopt":
		command => "tar -xzvf nlopt-2.2.4.tar.gz",
		cwd     => "${tools_dir}/treePL/deps",
		creates => "${tools_dir}/treePL/deps/nlopt-2.2.4",
		require => Exec["clone_treepl","install_openmpi"];
	"install_nlopt":
		command => "${tools_dir}/treePL/deps/nlopt-2.2.4/configure && make install",
		cwd		=> "${tools_dir}/treePL/deps/nlopt-2.2.4",
		creates => "/usr/local/include/nlopt.h",
		require => Exec["unzip_nlopt"];
	"compile_treepl":
		command => "${tools_dir}/treePL/src/configure && make",
		cwd     => "${tools_dir}/treePL/src",
		creates => "${tools_dir}/treePL/src/treePL",
		require => Exec["install_adolc","install_nlopt"];		
}
