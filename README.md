# Piscator 

Piscator is a [nexflow](https://www.nextflow.io/) computational pipeline for rapid evaluation primers or oligonucleotides used in the Polymerase Chain Reaction (PCR) for amplicon sequencing of microbial communities.

## **Getting Started**

Piscator is basically several utilities from the [primerprospector](http://pprospector.sourceforge.net) pipeline program, [EMBOSS](http://http://emboss.sourceforge.net/) utilities and several custom python scripts for connecting to a mysql data, analysing and plotting results. To get started clone/download this repository to your local machine and install the prerequisites in the list below. 


## Prerequisites
1.  Nextflow 
2.  primerprospector (python2.6)
3.  EMBOSS utilities
4.  MySQL

### Installation


1. First you must all the programs, tools, related dependencies required by the computer pipeline.  
   ```
   sudo apt-get install build-essential
   sudo apt update  && apt install -y \
       build-essential 
       cd-hit \
       clustalo \
       emboss \
       libfreetype6-dev \
       libpng-dev \
       libx11-dev \
       openjdk-8-jdk \
       python-pip \
       python-tk \
       python2.7
   ```

1. Install Nextflow more details from [nexflow](https://www.nextflow.io/)

1. Install the dependencies using using the file dependencies file
   `pip install  -r  main-requirements.txt`
   
1. Install **[Primerprospector](http://pprospector.sourceforge.net/install/install.html)**
   
   ```
     wget https://sourceforge.net/projects/pprospector/files/pprospector-1.0.1.tar.gz && \
     tar -zxvf  pprospector-1.0.1.tar.gz && \
     cd pprospector-1.0.1 && \
     pip install .
   ```
 
1. Install [RDP classifier](https://sourceforge.net/projects/rdp-classifier/)**
   * RDP classifier source files ``https://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.12.zip``
   * Add RDB to path ``echo "export PYTHONPATH=/home/pprospector/RDP/:$PYTHONPATH" >> /home/pprospector/.bashrc``
   * ``source /home/pprospector/.bashrc``

1. Python dependencies
   ```
   wget https://github.com/pycogent/pycogent/archive/1.5-release.tar.gz && \
       tar -zxvf   1.5-release.tar.gz && \
       cd pycogent-1.5-release && \
       pip install .
   ```
1. A workaround to address dependency issues. This is a temporary fix, Likely to change in the future.
   * Create an internal environment `virtualenv  python_virtualenv`
   * Activate the environment `source python_virtualenv/bin/activate`
   * Install the dependencies `pip install  -r ../env-requirements.txt`
   * Deactivate the environment `deactivate`

1. [MySQL](https://dev.mysql.com/downloads/mysql/) will need to installed and the hostname supplied to the `piscator.nf` and the password and username set in the `bin/init_primerDB.py` script


### **Running Piscator**
 *The main nextflow script file is `piscator.nf` and can be run using providing a tsv file containing the primers and changing the corresponding section in the main script.
 *The script may be run by invoking the command `nextflow  piscator.nf`
 *Alternatively the script may be called from anywhere by making it executable and adding the directory to the PATH environmental variable.


### **Running under Docker**
 * The pipeline has implemented to run on Docker; however, this is not fully tested. The [Dockerfile](https://github.com/PiscatorX/piscator-pipeline/blob/master/docker/Dockerfile) is provided.
