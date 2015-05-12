#!/bin/bash
#
# Setup the the box. This runs as root

apt-get -y update

apt-get -y install curl

dpkg -i puppetlabs-release-trusty.deb
sudo apt-get -y install puppet


# You can install anything you need here.
