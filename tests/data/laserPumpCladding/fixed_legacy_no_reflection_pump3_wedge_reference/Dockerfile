FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        git \
        libgl1 \
        libglu1-mesa \
        libhwloc-dev \
        libomp-dev \
        ninja-build \
        pkg-config \
        python3 \
        python3-dev \
        python3-pip \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/hase-venv \
    && /opt/hase-venv/bin/python -m pip install --no-cache-dir --upgrade pip setuptools wheel \
    && /opt/hase-venv/bin/python -m pip install --no-cache-dir \
        gmsh \
        matplotlib \
        numba \
        numpy \
        pybind11 \
        pytest \
        pyyaml \
        scikit-build-core \
        scipy \
        vtk

ENV PATH=/opt/hase-venv/bin:$PATH
WORKDIR /work
