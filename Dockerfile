FROM python:3.6.15-slim-buster AS base

# Download latest listing of available packages:
RUN apt-get -y update
# Upgrade already installed packages:
RUN apt-get -y upgrade &&  \
    apt-get update --fix-missing && \
    apt-get install -y git  build-essential 

FROM base AS python_build
RUN python -m pip install --upgrade pip && \
    pip install cython

FROM python_build AS whatshap
COPY . /opt
WORKDIR /opt
RUN pip install .
CMD ["/usr/local/bin/whatshap"]