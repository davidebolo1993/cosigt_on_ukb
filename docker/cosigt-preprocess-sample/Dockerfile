FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    ca-certificates \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    pigz \
    && rm -rf /var/lib/apt/lists/*

# Install samtools 1.22
RUN wget -q https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 && \
    tar -xf samtools-1.22.tar.bz2 && \
    cd samtools-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j $(nproc) && \
    make install && \
    cd .. && rm -rf samtools-1.22*

WORKDIR /work

# Copy preprocessing script into the image
COPY preprocess_sample_docker.sh /usr/local/bin/preprocess_sample_docker.sh
RUN chmod +x /usr/local/bin/preprocess_sample_docker.sh

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
