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

RUN git clone  https://github.com/PiscatorX/piscator-pipeline.git

COPY requirements.txt   .

RUN wget https://github.com/pycogent/pycogent/archive/1.5-release.tar.gz && \
    tar -zxvf   1.5-release.tar.gz && \
    cd pycogent-1.5-release && \
    python setup.py install

RUN wget https://sourceforge.net/projects/pprospector/files/pprospector-1.0.1.tar.gz && \
    tar -zxvf  pprospector-1.0.1.tar.gz && \
    cd pprospector-1.0.1 && \
    python setup.py install

RUN wget https://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.12.zip && \
    unzip rdp_classifier_2.12.zip 

ENV PYTHONPATH /docker/rdp_classifier_2.12:${PYTHONPATH}

ENV PATH ${PATH}:/docker/piscator-pipeline/bin

RUN pip install numpy==1.11.1 &&  pip install  -r  requirements.txt






