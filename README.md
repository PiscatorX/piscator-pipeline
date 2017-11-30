# Piscator 

Piscator is a [nexflow](https://www.nextflow.io/) computational pipeline for rapid evaluation primers or oligonucleotided used in the Polymerase Chain Reaction (PCR) for amplicon-based metagenomics.

## **Getting Started**

Piscator is basically several utilities from the [primerprospector](http://pprospector.sourceforge.net) pipeline program, [EMBOSS](http://http://emboss.sourceforge.net/) utilities and several custom python scripts for connecting to a mysql data, analysing and plotting results. To get started clone/download this repository to your local machine and install the prerequisites in the list below. 

### Prerequisites
1.  Nextflow 
2.  primerprospector
3.  EMBOSS utilities
4.  Mysql


###**Nextflow**

Nextflow requires Java to to run, for more details follow the details provided here
1. You will need to run the commands below to get it
* ``sudo apt-get update``
* ``sudo apt-get install default``
2. Now you can download next flow using the following command
````
This will generate ``nextflow`` file which you will have to move to your PATH. A safe place would be the /opt directory for example
``mv nextflow  /opt``
*Once you have moved the file you can test the installation by running ``nextflow`` in the terminal it should return a help message.  





<!-- ### Docker images, source files ### -->

<!-- 1. https://store.docker.com/community/images/nextflow/nextflow -->



### **Primerprospector**

1 .Firt you must install python dependensies.  Install  pythop-pip a utility for downloading python modules
``sudo apt-get install pip``


<!-- 1 .Numpy -->

<!-- 1. First you must install PyCogent -->


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


<!-- ###**EMBOSS utilities** -->

<!-- 1. EMBOSS package ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.5.7.tar.gz  -->


<!-- ### Installing -->






<!-- ``` -->
<!-- curl -s https://get.nextflow.io | bash -->
<!-- ``` -->
<!-- This will create a `nextflow` executable in your current directory.   -->

<!-- A step by step series of examples that tell you have to get a development env running -->

<!-- Say what the step will be -->

<!-- ``` -->
<!-- Give the example -->
<!-- ``` -->

<!-- And repeat -->

<!-- ``` -->
<!-- until finished -->
<!-- ``` -->

<!-- End with an example of getting some data out of the system or using it for a little demo -->

<!-- ## Running the tests -->

<!-- Explain how to run the automated tests for this system -->

<!-- ### Break down into end to end tests -->

<!-- Explain what these tests test and why -->

<!-- ``` -->
<!-- Give an example -->
<!-- ``` -->

<!-- ### And coding style tests -->

<!-- Explain what these tests test and why -->

<!-- ``` -->
<!-- Give an example -->
<!-- ``` -->

<!-- ## Deployment -->

<!-- Add additional notes about how to deploy this on a live system -->

<!-- ## Built With -->

<!-- * [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used -->
<!-- * [Maven](https://maven.apache.org/) - Dependency Management -->
<!-- * [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds -->

<!-- ## Contributing -->

<!-- Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us. -->

<!-- ## Versioning -->

<!-- We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).  -->

<!-- ## Authors -->

<!-- * **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth) -->

<!-- See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project. -->

<!-- ## License -->

<!-- This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details -->

## Acknowledgments

* PrimerProspector: de novo design and taxonomic analysis of PCR primers. William A. Walters (1,6), J. Gregory Caporaso (2,6), Christian L. Lauber (3), Donna Berg-Lyons (3), Noah Fierer (4), and Rob Knight (2,5). Bioinformatics (2011) 27(8): 1159-1161 first published online February 23, 2011 doi:http:///10.1093/bioinformatics/btr087

<!-- * Hat tip to anyone who's code was used -->
<!-- * Inspiration -->
<!-- * etc -->

