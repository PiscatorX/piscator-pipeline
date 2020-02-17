# Piscator 

Piscator is a [nexflow](https://www.nextflow.io/) computational pipeline for rapid evaluation primers or oligonucleotides used in the Polymerase Chain Reaction (PCR) for amplicon sequencing of microbial comunities.

## **Getting Started**

Piscator is basically several utilities from the [primerprospector](http://pprospector.sourceforge.net) pipeline program, [EMBOSS](http://http://emboss.sourceforge.net/) utilities and several custom python scripts for connecting to a mysql data, analysing and plotting results. To get started clone/download this repository to your local machine and install the prerequisites in the list below. 


## Prerequisites
1.  Nextflow 
2.  primerprospector (python2.6)
3.  EMBOSS utilities
4.  MySQL

### Nextflow

Nextflow requires Java to to run, for more details follow the details provided here
1. You will need to run the commands below to get it
 * ``sudo apt-get update``
 * ``sudo apt-get install default``
2. Now you can download next flow using the following command
 * ``wget -qO- get.nextflow.io | bash``
 * This will generate a ``nextflow`` file which you will have to move to your PATH.
 * Once you have moved the file you can test the installation by running ``nextflow`` in the terminal it should return a help message.  


### Primerprospector

1. First you must install python dependensies.  Install  pythop-pip a utility for downloading python modules
   * ``sudo apt-get install build-essential``
   *  ```sh
      apt update  && apt install -y \
      build-essential \
      clustalo \
      emboss \
      git \
      libfreetype6-dev \
      libpng-dev \
      libx11-dev \
      python-pip \
      python2.7 \
      python-tk \
      unzip \
      wget ```  
      
   * ``sudo apt-get install pip``

2. Install the dependencies using using the file dependencies file
   * ``pip install  -r  main-requirements.txt``
2. Install **[Primerprospector](http://pprospector.sourceforge.net/install/install.html)**
   
  ```shell
     wget https://sourceforge.net/projects/pprospector/files/pprospector-1.0.1.tar.gz && \
     tar -zxvf  pprospector-1.0.1.tar.gz && \
     cd pprospector-1.0.1 && \
     pip install . ```
 
5.  Install **[RDP classifier](https://sourceforge.net/projects/rdp-classifier/)**

 * RDP classifier source files ``https://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.12.zip``
 * ``echo "export PYTHONPATH=/home/pprospector/RDP/:$PYTHONPATH" >> /home/pprospector/.bashrc``
 * ``source /home/pprospector/.bashrc``

6. **Python dependencies**

 * PyCogent (ver. 1.5)  http://sourceforge.net/projects/pycogent/files/PyCogent/1.5/PyCogent-1.5.tgz/download (license: GPL)
 * Numpy (ver. 1.3.0)   http://sourceforge.net/projects/numpy/files/NumPy/1.3.0/numpy-1.3.0.tar.gz/download (license: BSD)
 * Matplotlib (ver. 0.98.5.3)  http://iweb.dl.sourceforge.net/project/matplotlib/OldFiles/matplotlib-0.98.5.3.tar.gz (license: BSD)




FROM  ubuntu:18.04

MAINTAINER Andrew Ndhlovu (drewxdvst@outlook.com)

ENV DEBIAN_FRONTEND=noninteractive

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
    python-tk \
    unzip \
    wget  

WORKDIR /docker

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

RUN pip install  -r  main-requirements.txt
ARG CACHEBUST=5
RUN git clone  https://github.com/PiscatorX/piscator-pipeline.git

ENV PATH /docker/piscator-pipeline/bin:\
/docker:\
${PATH}

RUN mkdir piscator-pipeline/python_virtualenv

RUN virtualenv  piscator-pipeline/python_virtualenv 

ENV VIRTUAL_ENV "/docker/piscator-pipeline/python_virtualenv"

ENV _OLD_VIRTUAL_PATH ${PATH}

ENV PATH ${VIRTUAL_ENV}/bin:${PATH}

RUN pip install  -r env-requirements.txt 

ENV PATH ${_OLD_VIRTUAL_PATH} 

RUN apt install -y cd-hit

ENV RDP_JAR_PATH /docker/rdp_classifier/rdp_classifier-2.0.jar

RUN apt install -y  openjdk-8-jdk
