# Based on https://hub.docker.com/r/litd/docker-cellranger
#FROM litd/docker-cellranger:v8.0.0

MAINTAINER Loréna Quatreville (lorena.quatreville@univ-amu.fr)

# Start from your existing Cell Ranger image
FROM litd/docker-cellranger:v8.0.1

# Set environment variables
ENV PATH=/opt/cellranger-8.0.1:$PATH


# Install dependencies for BCL Convert
RUN yum -y update && \
    yum -y install \
        curl \
        gcc \
        make \
        glibc \
        glibc-devel \
        zlib \
        zlib-devel \
        ca-certificates \
        rpm-build \
        gdb \
    && yum clean all

# Download BCL Convert v4.2.7 RPM from Illumina
RUN curl -L -o /tmp/bcl-convert-4.2.7.rpm \
    "https://s3.amazonaws.com/webdata.illumina.com/downloads/software/bcl-convert/bcl-convert-4.2.7-2.el8.x86_64.rpm" && \
    rpm -ivh /tmp/bcl-convert-4.2.7.rpm && \
    rm /tmp/bcl-convert-4.2.7.rpm

# Verify installation
RUN bcl-convert --version
RUN cellranger --version
