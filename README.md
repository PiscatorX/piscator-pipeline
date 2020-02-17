# Piscator 

Piscator is a [nexflow](https://www.nextflow.io/) computational pipeline for rapid evaluation primers or oligonucleotides used in the Polymerase Chain Reaction (PCR) for amplicon sequencing of microbial comunities.

## **Getting Started**

Piscator is basically several utilities from the [primerprospector](http://pprospector.sourceforge.net) pipeline program, [EMBOSS](http://http://emboss.sourceforge.net/) utilities and several custom python scripts for connecting to a mysql data, analysing and plotting results. To get started clone/download this repository to your local machine and install the prerequisites in the list below. 


## Prerequisites
1.  Nextflow 
2.  primerprospector (python2.6)
3.  EMBOSS utilities
4.  MySQL

### Installation

Nextflow requires Java to to run, more details from [nexflow](https://www.nextflow.io/)

1. First you must all the programs, tools, related dependensies required by the computer pipeline.  
   * ``sudo apt-get install build-essential``
   *  ```apt update  && apt install -y \
	     build-essential \
	     cd-hit \
	     clustalo \
	     emboss \
	     libfreetype6-dev \
	     libpng-dev \
	     libx11-dev \
	     openjdk-8-jdk \
	     python-pip \
	     python-tk \
	     python2.7 \
	     unzip \
	     wget ```  

2. Install the dependencies using using the file dependencies file
   * ```pip install  -r  main-requirements.txt```
   
3. Install **[Primerprospector](http://pprospector.sourceforge.net/install/install.html)**
   
  ```wget https://sourceforge.net/projects/pprospector/files/pprospector-1.0.1.tar.gz && \
     tar -zxvf  pprospector-1.0.1.tar.gz && \
     cd pprospector-1.0.1 && \
     pip install . ```
 
4.  Install **[RDP classifier](https://sourceforge.net/projects/rdp-classifier/)**

 * RDP classifier source files ``https://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.12.zip``
 * ``echo "export PYTHONPATH=/home/pprospector/RDP/:$PYTHONPATH" >> /home/pprospector/.bashrc``
 * ``source /home/pprospector/.bashrc``

5. **Python dependencies**

    pip install numpy==1.7.1
    wget https://github.com/pycogent/pycogent/archive/1.5-release.tar.gz && \
     	 tar -zxvf   1.5-release.tar.gz && \
     	 cd pycogent-1.5-release && \
     	 pip install .


 * PyCogent (ver. 1.5)  http://sourceforge.net/projects/pycogent/files/PyCogent/1.5/PyCogent-1.5.tgz/download (license: GPL)
 * Numpy (ver. 1.3.0)   http://sourceforge.net/projects/numpy/files/NumPy/1.3.0/numpy-1.3.0.tar.gz/download (license: BSD)
 * Matplotlib (ver. 0.98.5.3)  http://iweb.dl.sourceforge.net/project/matplotlib/OldFiles/matplotlib-0.98.5.3.tar.gz (license: BSD)
6. A workaround to address dependency issue. This is a temporary fix, Likely to change in the future.
   Create an internal environment ``virtualenv  python_virtualenv``
   Activate the environment ``source python_virtualenv/bin/activate``
   Install the dependencies  ``pip install  -r ../env-requirements.txt``
   
### **Running under Docker**
 * The pipeline has implemented to run on Docker; however, this is not fully tested. The [Dockerfile](https://github.com/PiscatorX/piscator-pipeline/blob/master/docker/Dockerfile) is provided.