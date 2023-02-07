# Get the base Ubuntu image from Docker Hub
FROM ubuntu:20.04

# Install GCC and dependencies
RUN apt-get -y update \
    && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata \
    && apt-get install -y git make cmake gcc g++

# Start building
COPY . /usr/src/rpfbwt
WORKDIR /usr/src/rpfbwt/build
RUN cmake -DENABLE_MIMALLOC=ON .. \
    && make -j

# Get binaries
WORKDIR /rpfbwt/bin
RUN cp \
    /usr/src/rpfbwt/build/rpfbwt32 \
    /usr/src/rpfbwt/build/rpfbwt64  \
    .
ENV PATH /rpfbwt/bin:$PATH