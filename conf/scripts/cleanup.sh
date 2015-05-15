#!/bin/bash

# Removing leftover leases and persistent rules
echo "cleaning up dhcp leases"
rm /var/lib/dhcp/*

# Make sure Udev doesn't block our network
echo "cleaning up udev rules"
if [ -d "/etc/udev/rules.d/70-persistent-net.rules" ]; then
	rm /etc/udev/rules.d/70-persistent-net.rules
	mkdir /etc/udev/rules.d/70-persistent-net.rules
fi
rm -rf /dev/.udev/
rm /lib/udev/rules.d/75-persistent-net-generator.rules

echo "Adding a 2 sec delay to the interface up, to make the dhclient happy"
echo "pre-up sleep 2" >> /etc/network/interfaces

# ADDED ADDITIONAL CLEANUP STEPS BELOW.
# SOURCE: http://vmassuchetto.github.io/2013/08/14/reducing-a-vagrant-box-size/

# Remove APT cache
apt-get clean -y
apt-get autoclean -y

# Remove APT files
find /var/lib/apt -type f | xargs rm -f

# Remove documentation files
if [ -d "/var/lib/doc" ]; then
	find /var/lib/doc -type f | xargs rm -f
fi

# Remove Virtualbox specific files
rm -rf /usr/src/vboxguest* /usr/src/virtualbox-ose-guest*

# Remove Linux headers
rm -rf /usr/src/linux-headers*

# Remove Unused locales (edit for your needs, this keeps only en* and pt_BR)
locales=`find /usr/share/locale/{af,am,ar,as,ast,az,be,bg,bn,bn_IN,br,bs,byn,ca,cs,csb,cy,da,de,dz,el,en_AU,en_CA,eo,es,et,eu,fa,fi,fo,fr,fur,ga,gez,gl,gu,haw,he,hi,hr,hu,hy,id,is,it,ja,ka,kk,km,kn,ko,kok,ku,ky,lt,lv,mg,mi,mk,ml,mn,mr,ms,mt,nb,ne,nl,nn,nso,oc,or,pa,pl,ps,ro,ru,rw,si,sk,sl,so,sq,sr,sr*latin,sv,sw,ta,te,th,ti,tig,tk,tl,tr,tt,ur,ve,vi,wa,wal,wo,xh,zh_HK,zh_CN,zh_TW,zu} -type d`
for loc in $locales; do
	rm -rf $loc
done

# Remove bash history
unset HISTFILE
rm -f /root/.bash_history
rm -f /home/vagrant/.bash_history
 
# Cleanup log files
find /var/log -type f | while read f; do echo -ne '' > $f; done;
