FROM ubuntu:18.04 as builder

ENV PACKAGES "g++ make cmake libboost-all-dev"
ENV SOURCE "/tmp/source"
ENV BUILD "/tmp/build"

RUN apt-get update &&\
    apt-get install -y ${PACKAGES}

COPY CMakeLists.txt freddi.ini ${SOURCE}/
COPY cpp ${SOURCE}/cpp

RUN mkdir -p ${BUILD} &&\
    cd ${BUILD} &&\
    cmake ${SOURCE} -DSTATIC_LINKING=TRUE -DNO_PYTHON_MODULE=TRUE &&\
    make install VERBOSE=1



FROM ubuntu:18.04

LABEL maintainer="Konstantin Malanchev <malanchev@sai.msu.ru>"

COPY --from=builder /usr/local/bin/freddi /bin/freddi
COPY --from=builder /usr/local/bin/freddi-ns /bin/freddi-ns

VOLUME /data
WORKDIR /
ENTRYPOINT ["/bin/freddi", "--dir=/data"]
