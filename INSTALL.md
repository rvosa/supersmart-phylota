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

Extended installation instructions
==================================

For installation instructions, please also refer to the 
[project website](http://www.supersmart-project.org/).
The SUPERSMART pipeline is designed to run on a virtual 
guest operating system within the 'VirtualBox' software, 
using the 'Vagrant' tool to insert the virtual guest operating
system into the VirtualBox.

This means that SUPERSMART is not actually installed on your own
operating system but runs on a virtual linux guest machine that is 
downloaded from a remote server onto your computer. Vagrant takes 
care of downloading the virtual image on which SUPERSMART runs. 
After downloading, not only SUPERSMART, but also all necessary 
dependencies (e.g. external software and libraries, databases) for 
the pipeline are automatically installed. It is thereby ensured 
that all software tools used in the pipeline run in a controlled
environment. 

You can download VirtualBox and Vagrant for all common operating 
systems at https://www.virtualbox.org/wiki/Downloads 
and http://www.vagrantup.com/downloads.html. You will find 
installation instructions on the respective websites.

Once Vagrant and VirtualBox are installed on your system, open a 
terminal ('command prompt' on Windows) and change to a suitable 
working directory. This directory will form the bridge between 
your operating system and the virtual guests (i.e. data files can
be exchanged between them), and some Vagrant-related files, including
a hidden folder, will be created here. It is therefore a good idea
to do this in a sensible subfolder, such as inside `supersmart-vagrant`
or something, rather than in your home folder, or desktop, or similar.
Once you are inside this folder, issue the following commands:

    $ vagrant init Naturalis/supersmart
    $ vagrant up

Vagrant will now download the virtual machine and SUPERSMART together 
with all necessary dependencies. Please not that this may take some time.
When the installation is finished, you can log into the virtual machine 
and change to the SUPERSMART working directory.

    $ vagrant ssh 
    $ cd $SUPERSMART_HOME

At this point it is always advisable to check whether there are updates 
to the pipeline source code. This is done by issuing the following 
command, which will download the latest code from the repository:

    $ git pull

For instructions on how to run the pipeline, please refer to the README 
file in the 'examples/' directory.

Logging into the guest machine on WINDOWS
=========================================

Please note that an ssh client is not installed by default on Windows 
systems. In order to log on to the virtual Linux guest machine, you 
need to download and run an ssh client, e.g. 'putty'; you will find
a download link at http://www.putty.org/ . Putty does not require any 
installation steps, you can just double click on the downloaded file 
and a dialog window will open asking you for the connection details. 
In the field 'Hostname (or IP address)', enter: `127.0.0.1`
In the field "Port", enter: `2222` and press 'Open'. A command prompt 
will open and ask you for a user name and password. For both, username 
and password, enter 'vagrant' to log in.
