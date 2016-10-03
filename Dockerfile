FROM ubuntu

MAINTAINER Konstantin Malanchev <malanchev@sai.msu.ru>

ENV PACKAGES "g++ make libboost-program-options-dev git"

RUN apt-get update && apt-get install -y ${PACKAGES}

WORKDIR /tmp
RUN git clone https://github.com/hombit/freddi.git
WORKDIR /tmp/freddi
RUN mkdir -p /usr/local/bin
RUN make install

RUN rm -r /tmp/freddi
RUN apt-get remove -y ${ENV} && apt-get autoremove

VOLUME /data
WORKDIR /
ENTRYPOINT ["/usr/local/bin/freddi", "--dir=/data"]
