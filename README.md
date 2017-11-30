# Piscator 

Piscator is a [nexflow](https://www.nextflow.io/) computational pipeline for rapid evaluation primers or oligonucleotided used in the Polymerase Chain Reaction (PCR) for amplicon-based metagenomics.

## **Getting Started**

Piscator is basically several utilities from the [primerprospector](http://pprospector.sourceforge.net) pipeline program, [EMBOSS](http://http://emboss.sourceforge.net/) utilities and several custom python scripts for connecting to a mysql data, analysing and plotting results. To get started clone/download this repository to your local machine and install the prerequisites in the list below. 

### Prerequisites
1.  Nextflow 
2.  primerprospector
3.  EMBOSS utilities
4.  Mysql


### Nextflow

Nextflow requires Java to to run, for more details follow the details provided here
1. You will need to run the commands below to get it
* ``sudo apt-get update``
* ``sudo apt-get install default``
2. Now you can download next flow using the following command
*``wget -qO- get.nextflow.io | bash``
This will generate a ``nextflow`` file which you will have to move to your PATH.
*Once you have moved the file you can test the installation by running ``nextflow`` in the terminal it should return a help message.  







### Primerprospector

1. Firt you must install python dependensies.  Install  pythop-pip a utility for downloading python modules
``sudo apt-get install pip``



1. **[Primerprospector](http://pprospector.sourceforge.net/install/install.html)**
   *Build-essential ``sudo apt-get install build-essential``
2.  Dependencies
   *Python 2.6 (src) (license: PSF)
   *PyCogent 1.5 (src) (license: GPL) 
   *Numpy 1.3.0 (src) (license: BSD)
   *Matplotlib 0.98.5.3 (src) (license: BSD)
   * The ``taxa_assignment_report.py`` script requires   
   
3. **Primerprospector installation**
* ``svn co https://pprospector.svn.sourceforge.net/svnroot/pprospector/trunk pprospector``
*  ``cd /home/pprospector``
*  ``tar -xvzf pprospector-1.0.1.tar.gz``
*  ``python setup.py install --install-scripts=/home/pprospector/bin/``
*  ``echo "export PATH=/home/pprospector/bin/:$PATH" >> /home/pprospector/.bashrc``
*  ``source /home/prospector/.bashrc``
*  ``cd /home/pprospector/tests/``

4.  **[RDP classifier](https://sourceforge.net/projects/rdp-classifier/) installation**
* RDP classifier source files ``https://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.12.zip``
* ``echo "export PYTHONPATH=/home/pprospector/RDP/:$PYTHONPATH" >> /home/pprospector/.bashrc``
* ``source /home/pprospector/.bashrc``

4. **Python dependencies**

*PyCogent (ver. 1.5)  http://sourceforge.net/projects/pycogent/files/PyCogent/1.5/PyCogent-1.5.tgz/download (license: GPL)
*Numpy (ver. 1.3.0)   http://sourceforge.net/projects/numpy/files/NumPy/1.3.0/numpy-1.3.0.tar.gz/download (license: BSD)
*Matplotlib (ver. 0.98.5.3)  http://iweb.dl.sourceforge.net/project/matplotlib/OldFiles/matplotlib-0.98.5.3.tar.gz (license: BSD)

