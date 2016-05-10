# This manifest installs the supersmart pipeline,
# for more information http://www.supersmart-project.org

# set user and default paths for storing data, tools and source code
$supersmart_home = "/supersmart"
$tools_dir       = "${supersmart_home}/tools"
$tools_bin_dir   = "${tools_dir}/bin"
$data_dir        = "${supersmart_home}/data"
$cran_url        = "http://cran.us.r-project.org"

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

# disable timeout for all provisioning operations
Exec { timeout => 0 }

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
	"libopenmpi-dev":  ensure => installed, require => Exec ["apt_update"];
	"openmpi-bin":     ensure => installed, require => Exec ["apt_update"];  
	"cpanminus":       ensure => installed, require => Exec ["apt_update"];  
	"phyml":           ensure => installed, require => Exec ["apt_update"];  
    "clang-3.4":       ensure => installed, require => Exec ["apt_update"];  
  }
  
  # create links for executables and data directories
  file {
	$supersmart_home:
	  ensure  => directory,
	  require => Exec["clone_supersmart"];
	$tools_dir:
	  ensure  => directory,
	  recurse => true;
	$tools_bin_dir:
	  ensure  => directory;
	$data_dir:
	  ensure  => directory;
    
	"muscle_link":
	  path    => "${tools_bin_dir}/muscle",
	  ensure  => link,
	  target  => "${tools_bin_dir}/muscle3.8.31_i86linux64",
	  require => Exec["unzip_muscle"];
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

	# install phylota database
	"dl_phylota_db":
	  command => "wget http://biovel.naturalis.nl/phylota.sqlite.gz",
	  cwd     => $data_dir,
	  creates => "${data_dir}/phylota.sqlite.gz",                        
	  require => [ File[ $data_dir ], Package[ 'wget' ], Exec[ 'clone_supersmart' ] ];    
	"unzip_phylota_db":		
	  command => "gunzip phylota.sqlite.gz",		
	  creates => "${data_dir}/phylota.sqlite",		
	  cwd     => $data_dir,		
	  require => Exec[ 'dl_phylota_db'];

	# install mafft
	"dl_mafft":
	  command   => "wget http://mafft.cbrc.jp/alignment/software/mafft-7.130-without-extensions-src.tgz",
	  cwd       => $tools_dir,
	  creates   => "${tools_dir}/mafft-7.130-without-extensions-src.tgz",
	  require   => [ File [ $tools_dir ] ];
	"unzip_mafft":
	  command   => "tar -xzvf mafft-7.130-without-extensions-src.tgz",
	  cwd       => $tools_dir,
	  creates   => "${tools_dir}/mafft-7.130-without-extensions",
	  require   => Exec[ 'dl_mafft' ];
	"make_nonroot_makefile_mafft":
	  command   => "sed -i -e 's|PREFIX = /usr/local|PREFIX = ${tools_dir}|g' Makefile",
	  cwd       => "${tools_dir}/mafft-7.130-without-extensions/core",
	  require   => Exec[ 'unzip_mafft' ];
	"install_mafft":
	  command   => "make && make install",
	  cwd       => "${tools_dir}/mafft-7.130-without-extensions/core",
	  creates   => "$tools_dir/bin/mafft",
	  require   => Exec[ 'make_nonroot_makefile_mafft' ];
    
	# install muscle multiple sequence alignment
	"download_muscle":
	  command   => "wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz",
	  cwd       => $tools_dir,
	  creates   => "${tools_dir}/muscle3.8.31_i86linux64.tar.gz",
	  require   => Package[ 'wget', 'tar' ];
	"unzip_muscle":
	  command   => "tar -xzvf muscle3.8.31_i86linux64.tar.gz -C ${tools_bin_dir}",
	  cwd       => $tools_dir,
	  creates   => "${tools_bin_dir}/muscle3.8.31_i86linux64",
	  require   => Exec["download_muscle"];
    
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
	  require => [ Exec["clone_examl"], Package["libopenmpi-dev","openmpi-bin"] ];
	"compile_parser":
	  command => "make -f Makefile.SSE3.gcc",
	  cwd     => "${tools_dir}/ExaML/parser",
	  creates => "${tools_dir}/ExaML/parser/parser",
	  require => [ Exec["clone_examl"], Package["libopenmpi-dev","openmpi-bin"] ];
    
	# install raxml
	"clone_raxml":
	  command => "git clone https://github.com/stamatak/standard-RAxML.git",
	  cwd     => $tools_dir,
	  creates => "${tools_dir}/standard-RAxML/",
	  require => Package[ 'git' ];
	"compile_raxml":
	  command => "make -f Makefile.SSE3.PTHREADS.gcc LD_LIBRARY_PATH=/usr/lib",
	  cwd     => "${tools_dir}/standard-RAxML/",
	  creates => "${tools_dir}/standard-RAxML/raxmlHPC-PTHREADS-SSE3",
	  require => Exec["clone_raxml"];
    
	# install exabayes
	"download_exabayes":
	  command => "wget http://sco.h-its.org/exelixis/material/exabayes/1.4.1/exabayes-1.4.1-linux-openmpi-sse.tar.gz",
	  cwd     => $tools_dir,
	  creates => "${tools_dir}/exabayes-1.4.1-linux-openmpi-sse.tar.gz",
	  require => Package["libopenmpi-dev","openmpi-bin"];
	"unzip_exabayes":
	  command => "tar -xvzf exabayes-1.4.1-linux-openmpi-sse.tar.gz",
	  cwd     => $tools_dir,
	  creates => "${tools_dir}/exabayes-1.4.1/",
	  require => Exec[ "download_exabayes" ];
	"configure_exabayes":		
	  command => "sh configure --enable-mpi",		
	  cwd     => "${tools_dir}/exabayes-1.4.1",		
      environment => ["CC=clang", "CXX=clang++"],
	  creates => "${tools_dir}/exabayes-1.4.1/exabayes/Makefile",		
	  require => [ Exec[ "unzip_exabayes" ], Package[ "clang-3.4" ] ];
	"compile_exabayes":		
	  command => "make",		
	  cwd     => "${tools_dir}/exabayes-1.4.1",		
	  creates => "${tools_dir}/exabayes-1.4.1/exabayes",		
	  require => Exec[ "configure_exabayes" ];
    
	# install supersmart
	"clone_supersmart":
	  command => "git clone https://github.com/naturalis/supersmart.git",
      cwd     => "/",
      creates => "/supersmart",
	  require => [ Package[ 'git' ] ];
        
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
	  require => [ Exec["clone_treepl"], Package["libopenmpi-dev","openmpi-bin"] ];
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
	  require => [ Exec["clone_treepl"], Package["libopenmpi-dev","openmpi-bin"] ];
	"install_nlopt":
	  command => "${tools_dir}/treePL/deps/nlopt-2.2.4/configure && make install",
	  cwd     => "${tools_dir}/treePL/deps/nlopt-2.2.4",
	  creates => "/usr/local/include/nlopt.h",
	  require => Exec["unzip_nlopt"];
	"compile_treepl":
	  command => "${tools_dir}/treePL/src/configure && make",
	  cwd     => "${tools_dir}/treePL/src",
	  creates => "${tools_dir}/treePL/src/treePL",
	  require => Exec["install_adolc","install_nlopt"];
	
	# install R packages
	"install_r_ape":
	  command => "echo 'install.packages(\"ape\",repos=\"${cran_url}\")' | R --vanilla",
	  require => Package['r-base', 'r-base-dev'];
	"install_r_phylosim":
	  command => "echo 'install.packages(\"phylosim\",repos=\"${cran_url}\")' | R --vanilla",
	  require => Exec['install_r_ape'];    	
	"install_r_phytools":
	  command => "echo 'install.packages(\"phytools\",repos=\"${cran_url}\")' | R --vanilla",
	  require => Exec['install_r_ape'];
	"install_r_phangorn":		
	  command => "echo 'install.packages(\"phangorn\",repos=\"${cran_url}\")' | R --vanilla",		
	  require => Exec['install_r_ape'];    

	  # install perl dependency libraries
	  "install_cpanm":
		command => "wget -O - https://cpanmin.us | sudo perl - App::cpanminus",
		require => Package['wget'];
	  "install_cpan_deps":
		command => "cpanm --notest --installdeps .",
		cwd     => "${supersmart_home}",
		require => Exec['install_cpanm','clone_supersmart'];
	  "install_bio_phylo":
		command => "cpanm --notest git://github.com/rvosa/bio-phylo.git",
		require => Exec['install_cpanm'];
	  "install_bioperl_live":
		command => "cpanm --notest git://github.com/bioperl/bioperl-live.git@v1.6.x",
		require => Exec['install_cpanm'];
	  "install_bioperl_run":
		command => "cpanm --notest git://github.com/bioperl/bioperl-run.git",
		require => Exec['install_cpanm'];    
  }
}
# Operations to shrink the size of the VM disk. Only run these cleanup operations when 
# provisioning a virtualbox image. 
class cleanup {
  if $virtual == 'virtualbox' {
	exec {
      
	  # clean up the VM
	  "clean_apt_get":
		command => "apt-get clean";
      "clean_meta":
		command => "rm MYMETA.*",	
		cwd     => "${supersmart_home}";
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
