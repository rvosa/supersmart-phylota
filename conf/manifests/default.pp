# == Class: supersmart
#
# This class installs the supersmart pipeline, 
# for more information http://www.supersmart-project.org
#
# === Parameters
#
# === Variables
#
# === Examples
#
# === Authors
#
# Author Name <rutgeraldo@gmail.com>
#
# === Copyright
#
# Copyright 2014 SUPERSMART project
#
#
class supersmart {

	# install the mysql server
	package { "mysql-server": ensure => installed }
	
	# install the mysql client
	package { "mysql": ensure => installed }	
	
	# set up the mysql daemon process
	service { "mysqld":
		enable  => true,
		ensure  => running,
		require => Package["mysql-server"],
	}

	# copy over the global mysql config file	
	file { "/var/lib/mysql/my.cnf":
		owner   => "mysql", 
		group   => "mysql",
		source  => "puppet:///mysql/my.cnf",
		notify  => Service["mysqld"],
		require => Package["mysql-server"],
	}
	
	# install the mysql config file
	file { "/etc/my.cnf":
		require => File["/var/lib/mysql/my.cnf"],
		ensure  => "/var/lib/mysql/my.cnf",
	}

	# set the mysql root password
	exec { "set-mysql-password":
		unless  => "mysqladmin -uroot -p$mysql_password status",
		path    => ["/bin", "/usr/bin"],
		command => "mysqladmin -uroot password $mysql_password",
		require => Service["mysqld"],
	}

	# this ensures all the perl dependencies are installed
	package { 'perl-doc': ensure => present }
	class { 'perl': require  => Package['perl-doc'] }
	perl::module { 'Bio::Phylo': require  => Package['perl-doc'] }
	perl::module { 'DBI': require  => Package['perl-doc'] }

	class { 'supersmart::instances': }

	vcsrepo { '/var/supersmart':
		ensure   => latest,
		provider => git,
		source   => 'https://github.com/naturalis/supersmart.git',
		require  => Class['supersmart::instances'],
	}
}