FROM  ubuntu:18.04

MAINTAINER Andrew Ndhlovu (drewxdvst@outlook.com)

RUN apt update  && apt install -y \
    build-essential \
    clustalo \
    emboss \
    git \
    python2.7 \
    python-pip \
    unzip \
    wget

WORKDIR /docker

RUN  git clone  https://github.com/PiscatorX/piscator-pipeline.git

ENV  PATH ${PATH}:/docker/piscator-pipeline/bin

COPY requirements.txt   .

RUN wget https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp-classifier%20v2.0/rdp_classifier_2.0.tar.gz

RUN tar -zxvf rdp_classifier_2.0.tar.gz 
   
ENV PYTHONPATH /docker/rdp_classifier

RUN pip install numpy==1.7.1

RUN wget https://github.com/pycogent/pycogent/archive/1.5-release.tar.gz && \
     tar -zxvf   1.5-release.tar.gz && \
     cd pycogent-1.5-release && \
     pip install .

RUN wget https://sourceforge.net/projects/pprospector/files/pprospector-1.0.1.tar.gz && \
     tar -zxvf  pprospector-1.0.1.tar.gz && \
     cd pprospector-1.0.1 && \
     pip install .

COPY *requirements.txt   /docker/

RUN apt update  && apt install -y \
    build-essential \
    clustalo \
    emboss \
    git \
    libfreetype6-dev \
    libpng-dev \
    libx11-dev \
    python-pip \
    python2.7 \
    unzip \
    wget  

RUN pip install  -r  main-requirements.txt

COPY test.py   /docker/

ENV PATH /docker:${PATH}




