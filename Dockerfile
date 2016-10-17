FROM ubuntu

MAINTAINER Konstantin Malanchev <malanchev@sai.msu.ru>

ENV PACKAGES "g++ make git libboost-program-options-dev"
ENV SOURCE "/tmp/freddi"

RUN apt-get update &&\
    apt-get install -y ${PACKAGES} &&\
    git clone https://github.com/hombit/freddi.git ${SOURCE} &&\
    cd ${SOURCE} &&\
    mkdir -p /usr/local/bin &&\
    make install LDLIBS=-static &&\
    rm -r ${SOURCE} &&\
    apt-get purge -y ${PACKAGES} &&\
    apt-get autoremove --purge -y

VOLUME /data
WORKDIR /
ENTRYPOINT ["/usr/local/bin/freddi", "--dir=/data"]
