# This file is read by docker to deploy the SUPERSMART environment
# on a virtual ubuntu host. Provisioning is done by puppet, using the
# 'default.pp' manifest.

# specify version of ubuntu
FROM ubuntu:14.04

# add scripts
ADD default.pp /
ADD intro.sh /

# install necessary packages
RUN apt-get -y update && apt-get install -y puppet nano

# provision virtual machine with SUPERSMART environment
RUN puppet apply /default.pp 

# set environment variables required by SUPERSMART
ENV PERL5LIB "$PERL5LIB:/supersmart/lib/"
ENV PATH "$PATH:/supersmart/script/"
ENV PATH "$PATH:/supersmart/tools/bin/"
ENV SUPERSMART_HOME "/supersmart/"
ENV EDITOR "nano"

# other environment variables
ENV SHELL "/bin/bash"

# set shared volume to share files with the virtual host
VOLUME [ "/shared" ]

# set directory on startup to shared dir
WORKDIR /shared

# automatically run intro script giving information on startup
RUN echo 'sh /intro.sh' >> /root/.bashrc
