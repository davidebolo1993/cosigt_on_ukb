FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install base dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    ca-certificates \
    git \
    pkg-config \
    cmake \
    curl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    pigz \
    parallel \
    r-base \
    python3-distutils \
    python3-dev \
    libjemalloc-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Go 1.19.1
RUN wget -q https://go.dev/dl/go1.19.1.linux-amd64.tar.gz && \
    tar -xf go1.19.1.linux-amd64.tar.gz -C /usr/local && \
    rm go1.19.1.linux-amd64.tar.gz
ENV PATH="/usr/local/go/bin:${PATH}"

# Install Rust and Cargo
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN rustup update
RUN cargo install --locked maturin

# Install samtools 1.22
RUN wget -q https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 && \
    tar -xf samtools-1.22.tar.bz2 && \
    cd samtools-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j $(nproc) && \
    make install && \
    cd .. && rm -rf samtools-1.22*

# Install BWA-MEM2 2.2.1
RUN wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    mv bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/ && \
    rm -rf bwa-mem2-2.2.1*

# Install bedtools 2.31.0
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static && \
    chmod +x bedtools.static && \
    mv bedtools.static /usr/local/bin/bedtools

# Install meryl
RUN wget -q https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz && \
    tar -xf meryl-1.4.1.Linux-amd64.tar.xz && \
    rm meryl-1.4.1.Linux-amd64.tar.xz
ENV PATH="/meryl-1.4.1/bin:${PATH}"

# Install kfilt (Go)
RUN git clone https://github.com/davidebolo1993/kfilt && \
    cd kfilt && \
    go mod init kfilt && \
    go mod tidy && \
    go build && \
    cp kfilt /usr/local/bin/ && \
    cd .. && rm -rf kfilt

# Install ODGI (CMake build) - use correct path
RUN git clone --recursive https://github.com/pangenome/odgi.git && \
    cd odgi && \
    cmake -H. -DCMAKE_BUILD_TYPE=Generic -Bbuild && \
    cmake --build build -- -j $(nproc) && \
    cp bin/odgi /usr/local/bin/odgi && \
    cd .. && rm -rf odgi

# Install gfainject (Rust)
RUN git clone https://github.com/AndreaGuarracino/gfainject && \
    cd gfainject && \
    cargo install --force --path . && \
    cd .. && rm -rf gfainject

# Install gafpack (Rust)
RUN git clone https://github.com/pangenome/gafpack.git && \
    cd gafpack && \
    cargo install --force --path . && \
    cd .. && rm -rf gafpack

# Install panplexity (Rust)
RUN git clone https://github.com/AndreaGuarracino/panplexity && \
    cd panplexity && \
    cargo install --force --path . && \
    cd .. && rm -rf panplexity

# Install cosigt (Go)
RUN git clone https://github.com/davidebolo1993/cosigt.git && \
    cd cosigt && \
    go mod init cosigt && \
    go mod tidy && \
    go build && \
    cp cosigt /usr/local/bin/ && \
    cd .. && rm -rf cosigt

# Install R packages (minimal set needed)
RUN R -e "install.packages(c('data.table', 'ggplot2', 'gridExtra', 'reshape2', 'rjson', 'dbscan', 'cluster'), repos='https://cloud.r-project.org/')"

# Copy all scripts and R files into the image
COPY process_region_docker.sh /usr/local/bin/process_region_docker.sh
COPY cluster.r /usr/local/bin/cluster.r
COPY coverage_outliers.r /usr/local/bin/coverage_outliers.r
RUN chmod +x /usr/local/bin/process_region_docker.sh
RUN chmod +x /usr/local/bin/*.r

WORKDIR /work

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

