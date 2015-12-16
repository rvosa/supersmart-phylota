FROM ubuntu:14.04

ADD default.pp /

RUN apt-get -y update && apt-get install -y puppet nano

RUN puppet apply /default.pp 

ENV PERL5LIB "$PERL5LIB:/supersmart/lib/"
ENV PATH "$PATH:/supersmart/script/"
ENV PATH "$PATH:/supersmart/tools/bin/"
ENV SUPERSMART_HOME "/supersmart/"
ENV EDITOR=nano

VOLUME [ "/shared" ]
WORKDIR /supersmart

RUN cd /supersmart
RUN /usr/bin/git pull
