# Set the base image to debian bionic
FROM ubuntu:18.04
USER root
ENV TZ=Asia/Tokyo
ENV CONDA_DIR $HOME/miniconda3
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# File Author / Maintainer
MAINTAINER Gergo Palfalvi <palfalvi.gergo@gmail.com>

RUN apt-get update && apt-get upgrade -y &&\
    apt-get install -y sudo \
    autoconf \
    automake \
    bamtools \
    bcftools \
    build-essential \
    cdbfasta \
    cmake \
    cpanminus \
    csh \
    default-jre \
    diamond-aligner \
    git \
    gnuplot \
    last-align \
    libatlas-base-dev \
    libbamtools-dev \
    libboost-all-dev \
    libdbd-mysql-perl \
    libfile-which-perl \
    libgsl-dev \
    libhts-dev \
    liblpsolve55-dev \
    libmysql++-dev \
    libparallel-forkmanager-perl \
    libpthread-stubs0-dev \
    libsqlite3-dev \
    libsuitesparse-dev \
    libyaml-perl \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    make \
    mc \
    mummer \
    python-matplotlib \
    python-numpy \
    python-pandas \
    python-pip \
    rsync \
    samtools \
    software-properties-common \
    time \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev




# Start working in /opt
WORKDIR /opt


# BRAKER installation
# HTSlib, braker, GeneMark-ES, Augustus
# GENEMARK-EX need a new link in every week or so
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
  URI::Escape \
  DB_File \
  threads &&\
  wget -c https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && tar -xf htslib-1.11.tar.bz2 && cd htslib-1.11 && ./configure && make && sudo make install && cd .. &&\
  wget -c https://github.com/Gaius-Augustus/BRAKER/archive/v2.1.5.tar.gz -O - | tar -xz &&\
  # wget -c https://github.com/lfaino/LoReAn/blob/master/third_party/software/gm_et_linux_64_4.48_3.60.tar.gz -O - | tar -xz &&\
  wget -c http://topaz.gatech.edu/GeneMark/tmp/GMtool_EcBpI/gmes_linux_64.tar.gz -O - | tar -xz &&\
  git clone https://github.com/Gaius-Augustus/Augustus && cd ./Augustus && sed -i 's/	cd bam2wig; make/	#cd bam2wig; make/g' ./auxprogs/Makefile && make clean && make && make install && cd .. &&\
  git clone https://github.com/splicebox/psiclass.git && cd psiclass && sudo make && cd .. &&\
  wget -c https://github.com/ruolin/strawberry/releases/download/v1.1.2/strawberry && chmod +x strawberry &&\
  git clone https://github.com/gpertea/stringtie && cd stringtie && make release && cd .. &&\
  wget -c https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.11.0/trinityrnaseq-v2.11.0.FULL.tar.gz -O - | tar -xz && cd trinityrnaseq-v2.11.0 && make && make install && cd .. &&\
  wget -c https://github.com/alexdobin/STAR/archive/2.7.7a.tar.gz -O - | tar -xz && cd STAR-2.7.7a/source && make STAR && cd ../.. &&\
  git clone https://github.com/lh3/minimap2 && cd minimap2 && make && cd .. &&\
  cpanm File::HomeDir

  #edta
  #portcullis
  #tama
  #unagi
  #taco
  #mikado
  #pasa

  #busco


ENV PATH=/opt/minimap2/:/opt/psiclass/:/opt/augustus-3.4.0/auxprogs/bam2hints:/opt/augustus-3.4.0/bin:/opt/augustus-3.4.0/scipts:/opt/canu/Linux-amd64/bin:/opt/BRAKER-2.1.5/scripts/:/opt/gmes_linux_64/:/opt/gmes_linux_64/ProtHint/bin/:$PATH
ENV AUGUSTUS_CONFIG_PATH=/opt/Augustus/config/
ENV AUGUSTUS_BIN_PATH=/opt/Augustus/bin/
ENV AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts/
