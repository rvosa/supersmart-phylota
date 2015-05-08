# This manifest installs the supersmart pipeline,
# for more information http://www.supersmart-project.org - the manifest is intended to be
# used in at least three contexts: i) to get the Travis continuous integration service up
# and running to do the unit tests; ii) to provision a vagrant box with all dependencies
# and shrink it as much as possible; iii) to install the pipeline natively on HPC 
# infrastructure such as an MPI compute cluster.

# set user and default paths for storing data, tools and source code
$username       = "vagrant"
$supersmart_dir = "/home/${username}/SUPERSMART"
$tools_dir      = "${supersmart_dir}/tools"
$tools_bin_dir  = "${tools_dir}/bin"
$src_dir        = "${supersmart_dir}/src"
$data_dir       = "${supersmart_dir}/data"
$cran_url       = "http://cran.us.r-project.org"

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

# This class contains the instructions for configuring and installing all dependencies,
# data, and the pipeline itself.
class install {

	# keep package information up to date
	exec {
		"apt_update":
		command => "/usr/bin/apt-get update"
	}

	# install packages.
	package {
		"mysql-server":    ensure => installed, require => Exec ["apt_update"];
		"mysql-client":    ensure => installed, require => Exec ["apt_update"];
		"ncbi-blast+":     ensure => installed, require => Exec ["apt_update"];
		"wget":            ensure => installed, require => Exec ["apt_update"];
		"tar":             ensure => installed, require => Exec ["apt_update"];
		"git":             ensure => installed, require => Exec ["apt_update"];
		"libarchive-dev":  ensure => installed, require => Exec ["apt_update"];
		"zlib1g-dev":      ensure => installed, require => Exec ["apt_update"];
		"autoconf":        ensure => installed, require => Exec ["apt_update"];
		"automake":        ensure => installed, require => Exec ["apt_update"];
		"libtool":         ensure => installed, require => Exec ["apt_update"];
		"build-essential": ensure => installed, require => Exec ["apt_update"];
		"curl":            ensure => installed, require => Exec ["apt_update"];
		"gzip":            ensure => installed, require => Exec ["apt_update"];
		"openjdk-6-jdk":   ensure => installed, require => Exec ["apt_update"];
		"subversion":      ensure => installed, require => Exec ["apt_update"];
		"r-base":          ensure => installed, require => Exec ["apt_update"];
		"r-base-dev":      ensure => installed, require => Exec ["apt_update"];

		# 2015-04-19 checking to see if we can get OpenMPI from package manager, as 
		# building from source takes so long that travis times out. Don't remember
		# why we weren't doing this in the first place. May have to revert this -- RAV
		"libopenmpi-dev":  ensure => installed, require => Exec ["apt_update"];
		"openmpi-bin":     ensure => installed, require => Exec ["apt_update"];  
	}


	# set up the mysql daemon process
	service {
		"mysql":
			enable  => true,
			ensure  => running,
			require => [ Package["mysql-server"], Exec["chown_phylota_db"] ];
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
			owner   => $username,
			recurse => true;
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
			target  => "${tools_dir}/exabayes-1.4.1/exabayes",
			require => Exec["compile_exabayes"];
		"parse_exabayes_link":
			path  => "${tools_bin_dir}/parse-exabayes",
			ensure  => link,
			target  => "${tools_dir}/exabayes-1.4.1/parser",
			require => Exec["compile_exabayes"];
		"exabayes_consense_link":
			path  => "${tools_bin_dir}/consense-exabayes",
			ensure  => link,
			target  => "${tools_dir}/exabayes-1.4.1/consense",
			require => Exec["compile_exabayes"];
		"raxml_link":
			path    => "${tools_bin_dir}/raxml",
			ensure  => link,
			target  => "${tools_dir}/standard-RAxML/raxmlHPC-PTHREADS-SSE3",
			require => Exec["compile_raxml"];
		"treepl_link":
			path    => "${tools_bin_dir}/treePL",
			ensure  => link,
			target  => "${tools_dir}/treePL/src/treePL",
			require => Exec["compile_treepl"];
	}

	# command line tasks
	exec {

		# set locale to US english to get rid of annoying perl warnings
		"set_locale_sh":
			command => "echo 'export LC_ALL=en_US.UTF-8' > set_locale.sh",
			cwd     => "/etc/profile.d",
			timeout => 0,        
			creates => "/etc/profile.d/set_locale.sh";

		# add bin directory for all required tools and smrt executables to PATH
		"make_bindir_sh":
			command => "echo 'export PATH=\$PATH:${tools_bin_dir}:${src_dir}/supersmart/script' > supersmart-tools-bin.sh",
			cwd     => "/etc/profile.d",
			timeout => 0,        
			creates => "/etc/profile.d/supersmart-tools-bin.sh";
		"make_bindir_csh":
			command => "echo 'setenv PATH \$PATH:${tools_bin_dir}:${src_dir}/supersmart/script' > supersmart-tools-bin.csh",
			cwd     => "/etc/profile.d",
			timeout => 0,        
			creates => "/etc/profile.d/supersmart-tools-bin.csh";

		# add default EDITOR environment variable
		"make_editor_sh":
			command => "echo 'export EDITOR=/etc/alternatives/editor' > supersmart-editor.sh",
			cwd     => "/etc/profile.d",
			timeout => 0,        
			creates => "/etc/profile.d/supersmart-editor.sh";
		"make_editor_csh":
			command => "echo 'setenv EDITOR=/etc/alternatives/editor' > supersmart-editor.csh",
			cwd     => "/etc/profile.d",
			timeout => 0,        
			creates => "/etc/profile.d/supersmart-editor.csh";

		# make phylota database
		"dl_phylota_dump":
			command => "wget http://biovel.naturalis.nl/phylota.tar.gz",
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
		"chown_phylota_db":
			command   => "chown -R -h mysql:mysql ${data_dir}/ && chmod 660 ${data_dir}/phylota/*",
			logoutput => true,
			timeout => 0,        
			require   => Exec[ 'unzip_phylota_dump' ];

		# install mafft
		"dl_mafft":
			command   => "wget http://mafft.cbrc.jp/alignment/software/mafft-7.130-without-extensions-src.tgz",
			cwd       => $tools_dir,
			creates   => "${tools_dir}/mafft-7.130-without-extensions-src.tgz",
			timeout => 0,        
			require   => [ File [ $tools_dir ] ];
		"unzip_mafft":
			command   => "tar -xzvf mafft-7.130-without-extensions-src.tgz",
			cwd       => $tools_dir,
			creates   => "${tools_dir}/mafft-7.130-without-extensions",
			timeout => 0,        
			require   => Exec[ 'dl_mafft' ];
		"make_nonroot_makefile_mafft":
			command   => "sed -i -e 's|PREFIX = /usr/local|PREFIX = ${tools_dir}|g' Makefile",
			cwd       => "${tools_dir}/mafft-7.130-without-extensions/core",
			timeout => 0,        
			require   => Exec[ 'unzip_mafft' ];
		"install_mafft":
			command   => "make && make install",
			cwd       => "${tools_dir}/mafft-7.130-without-extensions/core",
			creates   => "$tools_dir/bin/mafft",
			timeout => 0,        
			require   => Exec[ 'make_nonroot_makefile_mafft' ];

		# install muscle multiple sequence alignment
		"download_muscle":
			command   => "wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz",
			cwd       => $tools_dir,
			creates   => "${tools_dir}/muscle3.8.31_i86linux64.tar.gz",
			timeout => 0,        
			require   => Package[ 'wget', 'tar' ];
		"unzip_muscle":
			command   => "tar -xzvf muscle3.8.31_i86linux64.tar.gz -C ${tools_bin_dir}",
			cwd       => $tools_dir,
			creates   => "${tools_bin_dir}/muscle3.8.31_i86linux64",
			timeout => 0,        
			require   => Exec["download_muscle"];

		# 2015-04-19 checking to see if we can get OpenMPI from package manager, as 
		# building from source takes so long that travis times out. Don't remember
		# why we weren't doing this in the first place. May have to revert this -- RAV
	
		# 2015-05-06 So far this has not created any problems. We may NOT have to revert
		# this -- RAV

		# install openmpi
	#   "download_openmpi":
	#       command => "wget http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.5.tar.gz",
	#       cwd     => $tools_dir,
	#       creates => "${tools_dir}/openmpi-1.6.5.tar.gz",
	#       require => Package[ 'wget', 'tar' ];
	#   "unzip_openmpi":
	#       command => "tar -xvzf openmpi-1.6.5.tar.gz",
	#       cwd     => $tools_dir,
	#       creates => "${tools_dir}/openmpi-1.6.5",
	#       require => Exec['download_openmpi'];
	#   "install_openmpi":
	#       command => "${tools_dir}/openmpi-1.6.5/configure --prefix=/usr --disable-dlopen && make -j 4 install",
	#       creates => "/usr/bin/mpicc",
	#       cwd     => "${tools_dir}/openmpi-1.6.5",
	#       timeout => 0,
	#       require => Exec['unzip_openmpi'];

		# install phylip
		"download_phylip":
			command => "wget http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz",
			cwd     => $tools_dir,
			creates => "${tools_dir}/phylip-3.696.tar.gz",
			timeout => 0,        
			require => Package[ 'wget', 'tar' ];
		"unzip_phylip":
			command => "tar -xzvf ${tools_dir}/phylip-3.696.tar.gz",
			cwd     => $tools_dir,
			creates => "${tools_dir}/phylip-3.696/src/Makefile.unx",
			timeout => 0,        
			require => Exec["download_phylip"];
		"make_phylip":
			command => "make -f Makefile.unx consense",
			cwd     => "${tools_dir}/phylip-3.696/src/",
			creates => "${tools_dir}/phylip-3.696/src/consense",
			timeout => 0,        
			require => Exec["unzip_phylip"];

		# install examl & parser
		"clone_examl":
			command => "git clone https://github.com/stamatak/ExaML.git",
			cwd     => $tools_dir,
			creates => "${tools_dir}/ExaML/",
			timeout => 0,        
			require => Package[ 'git', 'zlib1g-dev' ];
		"compile_examl":
			command => "make -f Makefile.SSE3.gcc LD_LIBRARY_PATH=/usr/lib",
			cwd     => "${tools_dir}/ExaML/examl",
			creates => "${tools_dir}/ExaML/examl/examl",
			timeout => 0,        
			require => [ Exec["clone_examl"], Package["libopenmpi-dev","openmpi-bin"] ];
		"compile_parser":
			command => "make -f Makefile.SSE3.gcc",
			cwd     => "${tools_dir}/ExaML/parser",
			creates => "${tools_dir}/ExaML/parser/parser",
			timeout => 0,        
			require => [ Exec["clone_examl"], Package["libopenmpi-dev","openmpi-bin"] ];

		# install raxml
		"clone_raxml":
			command => "git clone https://github.com/stamatak/standard-RAxML.git",
			cwd     => $tools_dir,
			creates => "${tools_dir}/standard-RAxML/",
			timeout => 0,        
			require => Package[ 'git' ];
		"compile_raxml":
			command => "make -f Makefile.SSE3.PTHREADS.gcc LD_LIBRARY_PATH=/usr/lib",
			cwd     => "${tools_dir}/standard-RAxML/",
			creates => "${tools_dir}/standard-RAxML/raxmlHPC-PTHREADS-SSE3",
			timeout => 0,        
			require => Exec["clone_raxml"];

		# install exabayes
		"download_exabayes":
			command => "wget http://sco.h-its.org/exelixis/material/exabayes/1.4.1/exabayes-1.4.1-linux-openmpi-sse.tar.gz",
			cwd     => $tools_dir,
			creates => "${tools_dir}/exabayes-1.4.1-linux-openmpi-sse.tar.gz",
			timeout => 0,
			require => Package["libopenmpi-dev","openmpi-bin"];
		"unzip_exabayes":
			command => "tar -xvzf exabayes-1.4.1-linux-openmpi-sse.tar.gz",
			cwd     => $tools_dir,
			creates => "${tools_dir}/exabayes-1.4.1/",
			timeout => 0,        
			require => Exec[ "download_exabayes" ];
		"configure_exabayes":		
			command => "sh configure --enable-mpi",		
			cwd     => "${tools_dir}/exabayes-1.4.1",		
			creates => "${tools_dir}/exabayes-1.4.1/exabayes/Makefile",		
			timeout => 0,		
		require => Exec[ "unzip_exabayes" ];
		"compile_exabayes":		
			command => "make",		
			cwd     => "${tools_dir}/exabayes-1.4.1",		
			creates => "${tools_dir}/exabayes-1.4.1/exabayes",		
			timeout => 0,		
		require => Exec[ "configure_exabayes" ];

		# install supersmart
		"clone_supersmart":
			command => "git clone https://github.com/naturalis/supersmart.git",
			cwd     => $src_dir,
			creates => "${src_dir}/supersmart",
			user    => $username,
			timeout => 0,        
			require => [ File[ $src_dir ], Package[ 'git' ] ];
		"make_supersmart_sh":
			command => "echo 'export LD_LIBRARY_PATH=/usr/lib:/usr/lib64:/usr/local/lib' > supersmart.sh && echo 'export SUPERSMART_HOME=${src_dir}/supersmart' >> supersmart.sh && echo 'export PERL5LIB=\$PERL5LIB:\$SUPERSMART_HOME/lib' >> supersmart.sh",
			cwd     => "/etc/profile.d",
			timeout => 0,        
			creates => "/etc/profile.d/supersmart.sh";
		"make_supersmart_csh":
			command => "echo 'setenv LD_LIBRARY_PATH /usr/lib:/usr/lib64:/usr/local/lib' > supersmart.csh && echo 'setenv SUPERSMART_HOME ${src_dir}/supersmart' >> supersmart.csh && echo 'setenv PERL5LIB \$PERL5LIB:\$SUPERSMART_HOME/lib' >> supersmart.csh",
			cwd     => "/etc/profile.d",
			timeout => 0,        
			creates => "/etc/profile.d/supersmart.csh";

		# install BEAST
		"download_beast":
			command => "wget https://beast-mcmc.googlecode.com/files/BEASTv1.8.0.tgz",
			cwd     => $tools_dir,
			creates => "${tools_dir}/BEASTv1.8.0.tgz",
			timeout => 0,
			require => Package[ 'wget' ];
		"unzip_beast":
			command => "tar -xvzf BEASTv1.8.0.tgz",
			timeout => 0,        
			cwd     => "$tools_dir",
			creates => "${tools_dir}/BEASTv1.8.0/bin/beast",
			require => Exec[ 'download_beast' ];
		"symlink_beast":
			command => "ln -s ${tools_dir}/BEASTv1.8.0/bin/* .",
			timeout => 0,        
			cwd     => "/usr/local/bin/",
			creates => "/usr/local/bin/beast",
			require => Exec[ 'unzip_beast' ];

		# install beagle-lib
		"checkout_beagle_lib":
			command => "svn checkout http://beagle-lib.googlecode.com/svn/trunk/ beagle-lib",
			timeout => 0,        
			cwd     => $tools_dir,
			creates => "${tools_dir}/beagle-lib/autogen.sh",
			require => Package[ 'subversion', "openjdk-6-jdk" ];
		"generate_beagle_config":
			command => "sh autogen.sh",
			timeout => 0,        
			cwd     => "${tools_dir}/beagle-lib/",
			creates => "${tools_dir}/beagle-lib/configure",
			require => Exec[ 'checkout_beagle_lib' ];
		"build_beagle_lib":
			command => "sh configure && make install",
			timeout => 0,        
			cwd     => "${tools_dir}/beagle-lib/",
			creates => "/usr/local/lib/libhmsbeagle.so",
			require => Exec[ 'generate_beagle_config' ];

		# install treePL
		"clone_treepl":
			command => "git clone https://github.com/blackrim/treePL.git",
			timeout => 0,        
			cwd     => $tools_dir,
			creates => "${tools_dir}/treePL",
			require => Package[ 'git', 'autoconf', 'automake', 'libtool', 'build-essential', 'tar' ];
		"unzip_adolc":
			command => "tar -xvzf adol-c_git_saved.tar.gz",
			timeout => 0,        
			cwd     => "${tools_dir}/treePL/deps",
			creates => "${tools_dir}/treePL/deps/adol-c",
			require => [ Exec["clone_treepl"], Package["libopenmpi-dev","openmpi-bin"] ];
		"autoreconf_adolc":
			command => "autoreconf -fi",
			timeout => 0,        
			cwd     => "${tools_dir}/treePL/deps/adol-c",
			creates => "${tools_dir}/treePL/deps/adol-c/configure",
			require => Exec["unzip_adolc"];
		"install_adolc":
			command => "${tools_dir}/treePL/deps/adol-c/configure --prefix=/usr --with-openmp-flag=-fopenmp && make install",
			timeout => 0,        
			cwd     => "${tools_dir}/treePL/deps/adol-c",
			creates => "/usr/lib64/libadolc.so",
			require => Exec["autoreconf_adolc"];
		"unzip_nlopt":
			command => "tar -xzvf nlopt-2.2.4.tar.gz",
			timeout => 0,        
			cwd     => "${tools_dir}/treePL/deps",
			creates => "${tools_dir}/treePL/deps/nlopt-2.2.4",
			require => [ Exec["clone_treepl"], Package["libopenmpi-dev","openmpi-bin"] ];
		"install_nlopt":
			command => "${tools_dir}/treePL/deps/nlopt-2.2.4/configure && make install",
			timeout => 0,        
			cwd     => "${tools_dir}/treePL/deps/nlopt-2.2.4",
			creates => "/usr/local/include/nlopt.h",
			require => Exec["unzip_nlopt"];
		"compile_treepl":
			command => "${tools_dir}/treePL/src/configure && make",
			timeout => 0,        
			cwd     => "${tools_dir}/treePL/src",
			creates => "${tools_dir}/treePL/src/treePL",
			require => Exec["install_adolc","install_nlopt"];
		
		# install R packages
		"install_r_ape":
			command => "echo 'install.packages(\"ape\",repos=\"${cran_url}\")' | R --vanilla",
			timeout => 0,    	
			require => Package['r-base', 'r-base-dev'];
		"install_r_phangorn":
			command => "echo 'install.packages(\"phangorn\",repos=\"${cran_url}\")' | R --vanilla",
			timeout => 0,    	
			require => Exec['install_r_ape'];    	
		"install_r_phylosim":
			command => "echo 'install.packages(\"phylosim\",repos=\"${cran_url}\")' | R --vanilla",
			timeout => 0,    	
			require => Exec['install_r_ape'];    	
		"install_r_phytools":
			command => "echo 'install.packages(\"phytools\",repos=\"${cran_url}\")' | R --vanilla",
			timeout => 0,    	
			require => Exec['install_r_ape'];
	}

	# the Travis continuous integration system has special templates for testing different
	# perl versions in parallel. These templates use perlbrew and local::lib to hot swap 
	# different perl versions. For each version, cpanm can be used to install 
	# dependencies, and we do so using the commands in .travis.yml. For Travis we SHOULD 
	# NOT try to do that within this puppet manifest because then we end up installing 
	# dependencies system-wide, which means that the perlbrew that runs the unit tests 
	# cannot find them. Hence, we only install these here when FACTER variable $ci has 
	# not been set to 'travis'. To make sure that this variable IS set to 'travis', an
	# easy approach is to add the environment variable FACTER_ci=travis to the .travis.yml
	if $ci != "travis" {
		exec {
	
			# install perl dependency libraries
			"install_cpanm":
				command => "wget -O - https://cpanmin.us | sudo perl - App::cpanminus",
				timeout => 0,
				require => Package['wget'];
			"install_cpan_deps":
				command => "cpanm --notest --installdeps .",
				cwd     => "${src_dir}/supersmart",
				timeout => 0,
				require => Exec['install_cpanm','clone_supersmart'];
			"install_bio_phylo":
				command => "cpanm --notest git://github.com/rvosa/bio-phylo.git",
				timeout => 0,
				require => Exec['install_cpanm'];
			"install_bioperl_live":
				command => "cpanm --notest git://github.com/bioperl/bioperl-live.git@v1.6.x",
				timeout => 0,
				require => Exec['install_cpanm'];
			"install_bioperl_run":
				command => "cpanm --notest git://github.com/bioperl/bioperl-run.git",
				timeout => 0,
				require => Exec['install_cpanm'];
				
			# phylota database files need to be moved
			"mv_phylota_dump":
				command   => "mv ${data_dir}/phylota `mysql -s -N -uroot -p information_schema -e 'SELECT Variable_Value FROM GLOBAL_VARIABLES WHERE Variable_Name = \"datadir\"'`",
				logoutput => true,
				timeout => 0,        
				require   => Exec[ 'chown_phylota_db' ];                
			"grant_phylota_db":
				command   => "mysql -u root -e \"grant all on phylota.* to 'root'@'localhost';\"",
				logoutput => true,
				timeout => 0,        
				require   => Exec[ 'mv_phylota_dump' ]; 				
		}
	}
	else {	
		exec {
		
			# phylota database files need to be symlinked
			"mv_phylota_dump":
				command   => "ln -s ${data_dir}/phylota `mysql -s -N -uroot -p information_schema -e 'SELECT Variable_Value FROM GLOBAL_VARIABLES WHERE Variable_Name = \"datadir\"'`",
				logoutput => true,
				timeout => 0,        
				require   => Exec[ 'chown_phylota_db' ];                
			"grant_phylota_db":
				command   => "mysql -u root -e \"grant all on phylota.* to 'root'@'localhost';\"",
				logoutput => true,
				timeout => 0,        
				require   => Exec[ 'mv_phylota_dump' ];
		} 	
	}
}

# Operations to shrink the size of the VM disk. Only run these cleanup operations when 
# provisioning a virtualbox image. 
class cleanup {
	if $virtual == 'virtualbox' {
		exec {

			# clean up the VM
			"clean_apt_get":
				command => "apt-get clean",
				timeout => 0;		
		}
	}
}

# make sure all the cleanup happens last
stage { 'last': }
Stage['main'] -> Stage['last']
class { 'install':
      stage => main,
}
class { 'cleanup':
      stage => last,
}