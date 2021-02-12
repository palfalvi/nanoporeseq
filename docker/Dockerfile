# Set the base image to debian bionic
FROM ubuntu:18.04
USER root
ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# File Author / Maintainer
MAINTAINER Gergo Palfalvi <palfalvi.gergo@gmail.com>

RUN apt-get update && apt-get upgrade -y &&\
    apt-get install -y sudo \
    csh \
    cpanminus \
    git \
    make \
    build-essential \
    mummer \
    last-align \
    python-numpy \
    python-matplotlib \
    time \
    software-properties-common \
    gnuplot \
    zlib1g-dev \
    mc \
    wget \
    libatlas-base-dev \
    python-pip \
    python-pandas \
    cmake \
    default-jre \
    libbamtools-dev \
    bamtools \
    libmysql++-dev \
    libgsl-dev \
    libboost-all-dev \
    libsuitesparse-dev \
    liblpsolve55-dev \
    libsqlite3-dev \
    libboost-all-dev


# Start working in /opt
WORKDIR /opt


# BRAKER installation
# Perl modules
RUN cpanm \
  File::Spec::Functions \
  Hash::Merge \
  List::Util \
  Logger::Simple \
  Module::Load::Conditional \
  Parallel::ForkManager \
  POSIX \
  Scalar::Util::Numeric \
  YAML \
  Math::Utils \
  MCE::Mutex \
  Thread::Queue \
  threads &&\
  wget -c https://github.com/Gaius-Augustus/BRAKER/archive/v2.1.5.tar.gz -O - | tar -xz &&\
  wget -c http://topaz.gatech.edu/GeneMark/tmp/GMtool_BvMnP/gmes_linux_64.tar.gz -O - | tar -xz &&\
  git clone https://github.com/Gaius-Augustus/Augustus && cd Augustus && sudo make install &&\
  cd auxprogs/bam2hints && sudo make clean && cd /opt &&\
  cpanm File::HomeDir


ENV PATH=/opt/Augustus/auxprogs/bam2hints:/opt/Augustus/bin:/opt/Augustus/scipts:/opt/canu/Linux-amd64/bin:/opt/BRAKER-2.1.5/scripts/:/opt/gmes_linux_64/:/opt/gmes_linux_64/ProtHint/bin/:$PATH
ENV AUGUSTUS_CONFIG_PATH=/opt/Augustus/config/
ENV AUGUSTUS_BIN_PATH=/opt/Augustus/bin/
ENV AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts/
