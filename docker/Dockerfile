FROM ubuntu:xenial
MAINTAINER Hyeshik Chang <hyeshik@snu.ac.kr>

ENV LC_ALL C
ENV PATH /opt/tailseeker/bin:$PATH

RUN perl -pi -e 's,http://archive.ubuntu.com/ubuntu/,http://ftp.daum.net/ubuntu/,g' \
    /etc/apt/sources.list
RUN apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y install whiptail git pkg-config gcc wget make && \
    apt-get -y install libblas-dev liblapack-dev zlib1g-dev libbz2-dev gfortran && \
    apt-get -y install file python3 parallel python3-matplotlib python3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download tailseeker
RUN mkdir -p /opt; \
    cd /opt && \
    git clone https://github.com/hyeshik/tailseeker.git

# Install GSNAP/GMAP
RUN cd /tmp && \
    wget -q http://research-pub.gene.com/gmap/src/gmap-gsnap-2016-11-07.tar.gz && \
    tar -xzf gmap-gsnap-2016-11-07.tar.gz && \
    cd gmap-2016-11-07 && \
    env MAX_STACK_READLENGTH=512 ./configure --prefix=/opt/tailseeker && \
    make -j 8 && \
    make install && cd / && \
    rm -rf /tmp/gmap-gsnap-2016-11-07.tar.gz gmap-2016-11-07

# Install the STAR aligner
RUN cd /tmp && \
    wget -q https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz && \
    tar -xzf 2.5.2b.tar.gz && \
    cd STAR-2.5.2b/source && \
    make -j 8 STAR STARlong && \
    install -s -m 0755 -t /opt/tailseeker/bin STAR STARlong && \
    cd / && rm -rf /tmp/2.5.2b.tar.gz STAR-2.5.2b

# Install the AYB base-caller
RUN cd /tmp && \
    git clone https://github.com/hyeshik/AYB2.git && \
    cd AYB2/src && \
    make -j 8 && \
    install -s -m 0755 -t /opt/tailseeker/bin ../bin/AYB && \
    cd / && rm -rf /tmp/AYB2

# Install Python modules for tailseeker
RUN (wget -q https://bootstrap.pypa.io/ez_setup.py -O - | python3) && \
    easy_install pip && \
    rm -f setuptools*.zip && \
    pip install --no-cache-dir --upgrade cython && \
    pip install --no-cache-dir --upgrade snakemake colormath numpy scipy \
            pandas PyYAML biopython feather-format XlsxWriter

# Install htslib
RUN cd /tmp && \
    wget -q https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 && \
    tar -xjf htslib-1.3.1.tar.bz2 && \
    cd htslib-1.3.1 && \
    ./configure --prefix=/opt/tailseeker && make -j 8 && make install && \
    cd / && rm -rf /tmp/htslib-1.3.1.tar.bz2 /tmp/htslib-1.3.1

# Install samtools
RUN cd /tmp && \
    wget -q https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar -xjf samtools-1.3.1.tar.bz2 && \
    cd samtools-1.3.1 && \
    ./configure --prefix=/opt/tailseeker --with-htslib=/opt/tailseeker \
        --without-curses && \
    make -j 8 && make install && \
    cd / && rm -rf /tmp/samtools-1.3.1.tar.bz2 /tmp/samtools-1.3.1

# Install Heng Li's seqtk
RUN cd /tmp && \
    git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && make && \
    install -s -m 0755 -t /opt/tailseeker/bin seqtk && \
    cd / && rm -rf /tmp/seqtk

# Install bedtools
RUN cd /tmp && \
    wget -q https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz && \
    tar -xzf bedtools-2.26.0.tar.gz && \
    cd bedtools2 && \
    perl -pi -e 's/python/python3/' Makefile && \
    make -j 8 && \
    install -s -m 0755 -t /opt/tailseeker/bin bin/bedtools && \
    cd / && rm -rf /tmp/bedtools-2.26.0.tar.gz /tmp/bedtools2

# Install tailseeker 3
RUN cd /opt/tailseeker && \
    env TAILSEEKER_USE_AYB=yes TAILSEEKER_ANALYSIS_LEVEL=3 TAILSEEKER_USE_GSNAP=yes \
        TAILSEEKER_GSNAP_SENSITIVITY=yes TAILSEEKER_BINDIR=/opt/tailseeker/bin \
        PKG_CONFIG_PATH=/opt/tailseeker/lib/pkgconfig \
        sh setup.sh

# Prepare data directories.
ENV LD_LIBRARY_PATH /opt/tailseeker/lib
RUN mkdir -p -m 0777 /work /data
WORKDIR /work
ENTRYPOINT ["/opt/tailseeker/bin/tailseq-docker-wrap"]

# ex: ts=8 sts=4 sw=4 et
