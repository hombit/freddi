FROM ubuntu

MAINTAINER Konstantin Malanchev <malanchev@sai.msu.ru>

ENV COMPILETIME_PACKAGES "g++ make git"
ENV RUNTIME_PACKAGES "libboost-program-options-dev"
ENV SOURCE "/tmp/freddi"

RUN apt-get update &&\
    apt-get install -y ${COMPILETIME_PACKAGES} ${RUNTIME_PACKAGES} &&\
    git clone https://github.com/hombit/freddi.git ${SOURCE} &&\
    cd ${SOURCE} &&\
    mkdir -p /usr/local/bin &&\
    make install &&\
    rm -r ${SOURCE} &&\
    apt-get remove -y ${COMPILETIME_PACKAGES} &&\
    apt-get autoremove -y

VOLUME /data
WORKDIR /
ENTRYPOINT ["/usr/local/bin/freddi", "--dir=/data"]
