FROM shengqh/bioinfo:base

ENV BWA_VERSION="0.7.17"
RUN wget http://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2; \
    tar -jxvf bwa-${BWA_VERSION}.tar.bz2 ; \
    cd bwa-${BWA_VERSION}; \
    make; \
    cp bwa /usr/bin ; \
    cd .. ; \
    rm -rf bwa-${BWA_VERSION}

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir --upgrade parallel curl
    
ENV SVTOOL_VERSION="0.0.1"
RUN pip3 install git+git://github.com/zeissmania/sv_tool
RUN svtool -h
