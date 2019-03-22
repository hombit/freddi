#!/bin/bash

# Usage: .ci/install-boost.sh 69

BOOST_VERSION=1_${1}_0
BOOST_DOT_VERSION="1.${1}.0"

curl -LOJ https://dl.bintray.com/boostorg/release/${BOOST_DOT_VERSION}/source/boost_${BOOST_VERSION}.tar.gz
tar --gzip -xf boost_${BOOST_VERSION}.tar.gz
cd boost_${BOOST_VERSION}
./bootstrap.sh --with-libraries=program_options,python --with-python=python3
./b2 link=shared,static
./b2 install
