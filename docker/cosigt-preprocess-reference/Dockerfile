FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    ca-certificates \
    git \
    pigz \
    && rm -rf /var/lib/apt/lists/*

# Install meryl
RUN wget -q https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz && \
    tar -xf meryl-1.4.1.Linux-amd64.tar.xz && \
    rm meryl-1.4.1.Linux-amd64.tar.xz

ENV PATH="/meryl-1.4.1/bin:${PATH}"

WORKDIR /work

# Copy preprocessing script into the image
COPY preprocess_reference_docker.sh /usr/local/bin/preprocess_reference_docker.sh
RUN chmod +x /usr/local/bin/preprocess_reference_docker.sh

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
