FROM ubuntu:16.04

MAINTAINER Konstantin Malanchev <malanchev@sai.msu.ru>

ENV PACKAGES "g++ make libboost-all-dev"
ENV SOURCE "/tmp/freddi"

COPY ./ ${SOURCE}/

RUN apt-get update &&\
    apt-get install -y ${PACKAGES} &&\
    cd ${SOURCE} &&\
    mkdir -p /usr/local/bin &&\
    make install LDLIBS=-static &&\
    rm -r ${SOURCE} &&\
    apt-get purge -y ${PACKAGES} &&\
    apt-get autoremove --purge -y &&\
    apt-get clean -y &&\
    rm -rf /var/lib/apt/lists/* &&\
    truncate -s 0 /var/log/*log

VOLUME /data
WORKDIR /
ENTRYPOINT ["/usr/local/bin/freddi", "--dir=/data"]
