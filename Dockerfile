FROM ubuntu:18.04

# metadata
LABEL description="bito: Bayesian Inference of Trees via Optimization"

ENV APP_HOME /app
ENV DATA_HOME /data

# copy bito
WORKDIR $APP_HOME
COPY . $APP_HOME

# install dependencies
RUN apt-get update && \
  apt-get install -y \
  vim \
  git \
  wget \
  build-essential \
  gcc-7 \
  g++-7 \
  cmake \
  zlib1g-dev \
  python3-pip \
  python3-venv \
  python3-pytest \
  && apt-get clean

# install pip dependencies and build environment
RUN python3 -m venv /venv
RUN source /venv/bin/activate
RUN pip3 install numpy pandas
RUN echo "source /venv/bin/activate" >> ~/.bashrc

# install conda and build bito environment
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
#   /bin/bash miniconda.sh -b -p /opt/conda && \
#   rm miniconda.sh
# ENV PATH="/opt/conda/bin:${PATH}"
# RUN conda init
# WORKDIR $APP_HOME
# RUN conda env create -f environment.yml -n bito
# RUN echo "conda activate bito" >> ~/.bashrc
# RUN conda activate bito

# build bito
WORKDIR $APP_HOME
# symlink pytest
RUN if [ ! -e "pytest" ]; then ln -s $(which pytest-3) pytest; fi
RUN export PATH="$(pwd):$PATH"
RUN export CC=gcc-7
RUN export CXX=g++-7
RUN make clean
RUN make -j20 buildrelease

# set working directory
WORKDIR $DATA_HOME
