FROM ubuntu:14.04

ADD default.pp /

RUN apt-get -y update && apt-get install -y puppet

RUN puppet apply /default.pp 

ENV PERL5LIB "$PERL5LIB:/supersmart/lib/"
ENV PATH "$PATH:/supersmart/script/"
ENV SUPERSMART_HOME "/supersmart/"

VOLUME [ "/mydata" ]
WORKDIR /mydata

## CMD ~/supersmart/script/smrt

