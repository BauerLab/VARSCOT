# VARSCOT Dockerfile

FROM ubuntu:16.04

# General dependencies
RUN apt-get --yes update

RUN apt-get --yes install wget
RUN apt-get --yes install git
RUN apt-get --yes install gcc
RUN apt-get --yes install binutils
RUN apt-get --yes install make
RUN apt-get --yes install cmake

# Install Python and libraries
RUN apt-get --yes install python2.7
RUN apt-get install --yes python-setuptools
RUN apt-get install --yes python-pip
RUN apt-get install zlib1g-dev
RUN pip install sklearn>=0.19.0
RUN pip install numpy
RUN pip install scipy
RUN pip install pybedtools

# Install R
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -
RUN apt-get --yes update
RUN apt-get --yes install r-base r-base-dev

# Install required R libraries
RUN Rscript -e 'install.packages("randomForest", repos="http://R-Forge.R-project.org")'

# Lib dirs for Seqan and TUSCAN
RUN mkdir -p /app/lib

# Install Seqan
# https://github.com/seqan/seqan/tree/seqan-v2.3.2
WORKDIR /app/lib
RUN git clone --branch seqan-v2.4.0rc2 --depth 1 https://github.com/seqan/seqan.git

# Install TUSCAN
# https://github.com/BauerLab/TUSCAN/tree/master/TUSCAN%20model
WORKDIR /app/lib
RUN wget --no-verbose https://github.com/BauerLab/TUSCAN/archive/master.zip -O tuscan-master.zip
RUN unzip tuscan-master.zip
RUN rm tuscan-master.zip
RUN mv TUSCAN-master TUSCAN

# Build VARSCOT
WORKDIR /app
COPY . /app/

RUN mkdir -p /app/build/read_mapping_build
WORKDIR /app/build/read_mapping_build
RUN cmake ../../read_mapping -DCMAKE_BUILD_TYPE=Release
RUN make

RUN mkdir -p /app/build/variant_processing_build
WORKDIR /app/build/variant_processing_build
RUN cmake ../../variant_processing -DCMAKE_BUILD_TYPE=Release
RUN make

# Prepare to make use of volumes in pipeline command 
RUN mkdir -p /vcf
VOLUME /data/vcf
RUN mkdir -p /bed
VOLUME /data/bed
RUN mkdir -p /genome
VOLUME /data/genome
RUN mkdir -p /index
VOLUME /data/index
RUN mkdir -p /output
VOLUME /output
RUN mkdir -p /tmpdir
VOLUME /tmpdir

# Set pipeline entry point script
WORKDIR /app
COPY ./entrypoint.sh .
ENTRYPOINT ["/app/entrypoint.sh"]
