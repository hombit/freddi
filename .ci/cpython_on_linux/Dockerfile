ARG PY_VERSION
FROM python:${PY_VERSION}

ENV PACKAGES "g++ make curl"
ENV PYPKG "numpy scipy cmake"
ENV SOURCE "/tmp/source"
ENV BOOST_VERSION "69"

RUN apt-get update &&\
    apt-get install -y ${PACKAGES}

RUN pip3 install ${PYPKG}

COPY .ci/install-boost.sh /install-boost.sh

RUN /install-boost.sh ${BOOST_VERSION} &&\
    ldconfig

COPY . ${SOURCE}/

WORKDIR ${SOURCE}

RUN python setup.py install
