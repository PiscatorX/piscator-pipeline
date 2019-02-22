FROM  ubuntu:18.04

MAINTAINER Andrew Ndhlovu (drewxdvst@outlook.com)

RUN apt update  && apt install -y \
    git \
    python2.7 \
    python-pip \
    wget 

WORKDIR /docker

RUN git clone  https://github.com/PiscatorX/piscator-pipeline.git

COPY requirements.txt   .

RUN pip install numpy==1.11.1 &&  pip install  -r  requirements.txt


    

