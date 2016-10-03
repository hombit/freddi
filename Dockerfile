FROM ubuntu

MAINTAINER Konstantin Malanchev <malanchev@sai.msu.ru>

RUN apt-get update && apt-get install -y g++ make libboost-program-options-dev git

WORKDIR /tmp
RUN git clone https://github.com/hombit/freddi.git
WORKDIR /tmp/freddi
RUN mkdir -p /usr/local/bin
RUN make install
RUN rm -r /tmp/freddi

VOLUME /data
WORKDIR /
ENTRYPOINT ["/usr/local/bin/freddi"]
