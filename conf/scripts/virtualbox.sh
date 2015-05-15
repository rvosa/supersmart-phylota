#!/bin/bash

# Bail if we are not running inside VirtualBox.
if [[ `facter virtual` != "virtualbox" ]]; then
    exit 0
fi

# this is ubuntu-specific, which we are using so we can get away with
# doing it like this. The commands below mount the generic VBoxGuest 
# ISO and run a great big shell script that tries to compile the guest
# utils from source. Unfortunately that fails miserably. Hopefully this
# doesn't.
apt-get -y install virtualbox-guest-utils

# mkdir -p /mnt/virtualbox
# mount -o loop /home/vagrant/VBoxGuest*.iso /mnt/virtualbox
# sh /mnt/virtualbox/VBoxLinuxAdditions.run
# ln -s /opt/VBoxGuestAdditions-*/lib/VBoxGuestAdditions /usr/lib/VBoxGuestAdditions
# umount /mnt/virtualbox
# rm -rf /home/vagrant/VBoxGuest*.iso
