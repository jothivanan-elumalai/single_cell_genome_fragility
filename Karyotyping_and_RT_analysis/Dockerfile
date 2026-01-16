# Use Ubuntu 24.04 as base
FROM ubuntu:24.04

LABEL authors="Jothivanan Elumalai <jothivanan.elumalai[at]riken.jp and Hisashi Miura <hisashi.miura[at]riken.jp>" \
      description="Docker image for replication timing and karyotype analysis at the single-cell level"            

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    cmake \
    make \
    libhdf5-dev \
    language-pack-en \
    autoconf \
    automake \
    libtool \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev

# Set locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

# Install Miniforge3 v25.3.1
RUN wget https://github.com/conda-forge/miniforge/releases/download/25.3.1-0/Miniforge3-25.3.1-0-Linux-x86_64.sh && bash Miniforge3-25.3.1-0-Linux-x86_64.sh -p /miniforge3 -b
ENV PATH=/miniforge3/bin:$PATH
RUN rm Miniforge3-25.3.1-0-Linux-x86_64.sh

# Add necessary channels and set channel priority
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda config --set channel_priority strict

# create screpliseq conda environment using the provided YAML file
# Ensure the screpliseq.yml file is in the same directory as this Dockerfile
WORKDIR /usr/local/bin
COPY screpliseq.yml .
RUN conda env create -n screpliseq -f screpliseq.yml -y

# Activate screpliseq environment by default
SHELL ["/bin/bash", "-c"]
RUN echo "source activate screpliseq" >> ~/.bashrc
COPY entrypoint.sh .
RUN chmod +x entrypoint.sh
ENTRYPOINT ["./entrypoint.sh"]

# Install samstat v2.2.3
RUN git clone https://github.com/TimoLassmann/samstat.git && cd samstat && mkdir build && cd build && cmake .. && make
ENV PATH=/usr/local/bin/samstat/build/src/:$PATH

# Copy scripts and set permissions
RUN mkdir -p scripts
COPY scripts/*.sh scripts/
RUN mkdir -p util
COPY util/ util/
ENV PATH=/usr/local/bin/scripts/:$PATH
ENV PATH=/usr/local/bin/util/:$PATH
RUN chmod -R +x scripts
RUN chmod -R +x util
CMD ["pipeline-list.sh"]