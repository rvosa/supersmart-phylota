# This manifest installs the supersmart pipeline, 
# for more information http://www.supersmart-project.org

# set up the mysql daemon process
# 	service { "mysqld":
# 		enable  => true,
# 		ensure  => running,
# 		require => Package["mysql-server"],
# 	}
# 
# copy over the global mysql config file	
# 	file { "/var/lib/mysql/my.cnf":
# 		owner   => "mysql", 
# 		group   => "mysql",
# 		source  => "puppet:///mysql/my.cnf",
# 		notify  => Service["mysqld"],
# 		require => Package["mysql-server"],
# 	}
# 	
# install the mysql config file
# 	file { "/etc/my.cnf":
# 		require => File["/var/lib/mysql/my.cnf"],
# 		ensure  => "/var/lib/mysql/my.cnf",
# 	}
# 
# set the mysql root password
# 	exec { "set-mysql-password":
# 		unless  => "mysqladmin -uroot -p$mysql_password status",
# 		path    => ["/bin", "/usr/bin"],
# 		command => "mysqladmin -uroot password $mysql_password",
# 		require => Service["mysqld"],
# 	}
# vcsrepo { 
# 	"/usr/local/src/supersmart":
# 		ensure   => latest,
# 		provider => git,
# 		source   => 'https://github.com/naturalis/supersmart.git';
# 
# }

Exec {
	path => [ 
		"/usr/local/sbin", 
		"/usr/local/bin",
		"/usr/sbin",
		"/usr/bin",
		"/sbin",
		"/bin",
		"/usr/lib64/openmpi/bin/",
	]
}

package {
	"mysql-server":             ensure => installed;
	"mysql":                    ensure => installed;
	"ncbi-blast+.x86_64":       ensure => installed;
	"wget":                     ensure => installed;
	"tar":                      ensure => installed;
	"perl-DBI":                 ensure => installed;
	"perl-DBD-MySQL":           ensure => installed;
	"perl-DBIx-Class":          ensure => installed;
	"perl-Template-Toolkit":    ensure => installed;
	"perl-JSON":                ensure => installed;
	"perl-Moose":               ensure => installed;
	"perl-XML-Twig":            ensure => installed;
	"perl-HTML-Parser":         ensure => installed;
	"perl-Config-Tiny":         ensure => installed;
	"perl-bioperl":             ensure => installed;
	"openmpi":                  ensure => installed;
	"git":                      ensure => installed;
	"openmpi-devel":            ensure => installed;
	"zlib-devel":               ensure => installed;
	"autoconf":                 ensure => installed;
	"automake":                 ensure => installed;
	"libtool":                  ensure => installed;
	"gcc-c++":                  ensure => installed;	
}

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
}

exec {

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
		
	# install perl package Parallel::MPI::Simple
	"download_parallel_mpi_simple":
		command => "wget http://search.cpan.org/CPAN/authors/id/A/AJ/AJGOUGH/Parallel-MPI-Simple-0.10.tar.gz",
		cwd     => "/usr/local/src",
		creates => "/usr/local/src/Parallel-MPI-Simple-0.10.tar.gz",
		require => Package[ 'wget', 'tar', 'openmpi', 'openmpi-devel' ];		
	"unzip_parallel_mpi_simple":
		command => "tar -xzvf /usr/local/src/Parallel-MPI-Simple-0.10.tar.gz",
		cwd     => "/usr/local/src",		
		creates => "/usr/local/src/Parallel-MPI-Simple-0.10/Makefile.PL",
		require => Exec["download_parallel_mpi_simple"];
	"make_makefile_parallel_mpi_simple":
		command => "perl Makefile.PL CCFLAGS='-I/usr/include/openmpi-x86_64/'",
		cwd     => "/usr/local/src/Parallel-MPI-Simple-0.10/",
		creates => "/usr/local/src/Parallel-MPI-Simple-0.10/Makefile",
		require => Exec["unzip_parallel_mpi_simple"];
	"make_install_parallel_mpi_simple":
		command => "make install",
		cwd     => "/usr/local/src/Parallel-MPI-Simple-0.10/",		
		creates => "/usr/local/src/Parallel-MPI-Simple-0.10/Simple.so",
		require => Exec["make_makefile_parallel_mpi_simple"];		
	
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
		require => Package[ 'git', 'openmpi-devel', 'zlib-devel' ];
	"compile_examl":
		command => "make -f Makefile.SSE3.gcc",
		cwd     => "/usr/local/src/ExaML/examl",
		creates => "/usr/local/src/ExaML/examl/examl",
		require => Exec["clone_examl"];
	"compile_parser":
		command => "make -f Makefile.SSE3.gcc",
		cwd     => "/usr/local/src/ExaML/parser",
		creates => "/usr/local/src/ExaML/parser/parser",
		require => Exec["clone_examl"];
	
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
		require => Exec["clone_treepl"];
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
		require => Exec["clone_treepl"];
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