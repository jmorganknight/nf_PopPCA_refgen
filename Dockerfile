FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install required runtime tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    bzip2 \
    ca-certificates \
    bcftools \
    tabix \
    python3 \
    gawk \
    && rm -rf /var/lib/apt/lists/*

# Install PLINK 1.9
RUN wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip -O /tmp/plink1.zip \
    && unzip -q /tmp/plink1.zip -d /tmp/plink1 \
    && mv /tmp/plink1/plink /usr/local/bin/plink \
    && rm -rf /tmp/plink1.zip /tmp/plink1 \
    && chmod +x /usr/local/bin/plink

# Install PLINK 2.0
RUN wget -q https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip -O /tmp/plink2.zip \
    && unzip -q /tmp/plink2.zip -d /tmp/plink2 \
    && mv /tmp/plink2/plink2 /usr/local/bin/plink2 \
    && rm -rf /tmp/plink2.zip /tmp/plink2 \
    && chmod +x /usr/local/bin/plink2

WORKDIR /workspace

CMD ["/bin/bash"]
