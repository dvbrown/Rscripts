# Idenitfying deletions from mitochondrial DNA



## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisities

You will need:

A python installation. This pipeline has been designed with python 2.7. 
The python binary is located in

```
/cm/shared/apps/python/2.7.3/bin/python2.7
```
You will need to install the ruffus python module. To do this you need to install it in your local folder.
```
mkdir /home/dbrown0/local/lib/python2.7/site-packages
python2.7 -m easy_install --prefix=$HOME/local ruffus --upgrade

```
When you launch python2.7 ruffus should be installed

Use conda and the bioinformatics environment I made

source activate bioinformatics

```

### Installing

The pipeline code is available from github. Copy and paste this command if you have github command line tool installed: 
```
git clone https://github.com/dvbrown/Pipelines/tree/master/Ruffus_python/160306_mtDNApipelines
```
Alternatively you can use the github website and clone to your computer using a graphical user interface.

## Configuring global parameters for the pipeline
The paths to various software need to be set at the very top of the pindel_commands.py file.

```
Path to reference genome
```
```
Path to folder containing pindel executables
```
```
Path to GATK jar file
```

## Preparing files prior to running the pipeline


## Authors

* ** Insert authors here ** - *Initial work* - 
## License

This project is licensed under the GNU GPL License.

## Acknowledgments
