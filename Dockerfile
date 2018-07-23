FROM ubuntu:18.04 as builder

ENV PACKAGES "g++ make cmake libboost-all-dev"
ENV SOURCE "/tmp/source"
ENV BUILD "/tmp/build"

RUN apt-get update &&\
    apt-get install -y ${PACKAGES}

COPY ./ ${SOURCE}/

RUN mkdir -p ${BUILD} &&\
    cd ${BUILD} &&\
    cmake ${SOURCE} -DSTATIC_LINKING=TRUE &&\
    make install


FROM alpine

LABEL maintainer="Konstantin Malanchev <malanchev@sai.msu.ru>"

COPY --from=builder /usr/local/bin/freddi /bin/freddi

VOLUME /data
WORKDIR /
ENTRYPOINT ["/bin/freddi", "--dir=/data"]
