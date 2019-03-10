FROM ubuntu:18.04 as builder

ENV PACKAGES "g++ make cmake libboost-all-dev"
ENV SOURCE "/tmp/source"
ENV BUILD "/tmp/build"

RUN apt-get update &&\
    apt-get install -y ${PACKAGES}

RUN mkdir -p ${BUILD}
WORKDIR ${BUILD}

COPY ./ ${SOURCE}/

RUN cmake ${SOURCE} -DSTATIC_LINKING=TRUE &&\
    make install -j

RUN make test


FROM alpine

LABEL maintainer="Konstantin Malanchev <malanchev@sai.msu.ru>"

COPY --from=builder /usr/local/bin/freddi /bin/freddi
COPY --from=builder /usr/local/bin/freddi-ns /bin/freddi-ns

VOLUME /data
WORKDIR /
ENTRYPOINT ["/bin/freddi", "--dir=/data"]
