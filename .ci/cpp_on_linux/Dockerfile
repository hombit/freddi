FROM ubuntu:16.04

ENV PACKAGES "g++ make cmake curl"
ENV SOURCE "/tmp/source"
ENV BUILD "/tmp/build"
ENV BOOST_VERSION "69"

RUN apt-get update &&\
    apt-get install -y ${PACKAGES}

COPY . ${SOURCE}/

RUN ${SOURCE}/.ci/install-boost.sh ${BOOST_VERSION} &&\
    ldconfig

RUN mkdir -p ${BUILD}
WORKDIR ${BUILD}

RUN cmake ${SOURCE} &&\
    make VERBOSE=1 &&\
    make install VERBOSE=1
