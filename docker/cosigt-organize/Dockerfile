FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y git && \
    rm -rf /var/lib/apt/lists/*

# Install PyYAML
RUN pip install --no-cache-dir pyyaml

WORKDIR /work

# Copy organize.py into the image
COPY organize.py /usr/local/bin/organize.py
RUN chmod +x /usr/local/bin/organize.py

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
