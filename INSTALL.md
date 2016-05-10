Requirements 
============

- an installation of 'VirtualBox'
- an installation of 'Vagrant'
- a connection to the internet
- an SSH client

Quick installation guide (see below for extended version) 
=========================================================

- If not already installed, download and 
   install [VirtualBox](https://www.virtualbox.org/wiki/Downloads)
- If not already installed, download and
   install [Vagrant](http://www.vagrantup.com/downloads.html).
- Change to a suitable working directory and issue the following commands:
- `vagrant init Naturalis/supersmart`
- `vagrant up`
- Log in to the virtual guest machine using `vagrant ssh`
- SUPERSMART is installed in the directory `/home/vagrant/SUPERSMART`, i.e. $SUPERSMART_HOME
- The working directory where you typed `vagrant up` on the host is `/vagrant` inside the VM

Extended installation instructions
==================================
Extended installation instructions can be found on our [wiki page](https://github.com/naturalis/supersmart/wiki/Installation-instructions).
