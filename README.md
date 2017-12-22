# Regulatory Annotations Comparisons
This repository contains code to reproduce results of the manuscript:

Each directory contains code for analyses specific to different figures of the paper. Everything is run using [snakemake](http://snakemake.readthedocs.io/en/latest/. These directories follow this general pattern:
```	
├── .gitignore
├── README.md
├── config.yml : Cluster config specifications 
├── scripts : Scripts for analyses
│   ├── script1.py
│   └── script2.R
├── envs : Environment specification 
│   └── myenv.yml
├── run.sh : Commands for dry run, executing Snakemake and creating workflow DAGs
└── Snakefile : To run the workflow
```

# Requirements

The analyses use the following software:

* Python 3.6.3
* Snakemake 4.2.0
* R 3.3.3	
* bedtools v2.26.0
* GREGOR 1.2.1
* GAT

To setup these pre-requisites, use the Anaconda/Miniconda Python3 distribution. The Conda package manager is used to obtain and deploy the defined software packages in the specified versions. These instructions are for the Linux platform
	
# Step 1: Install [Anaconda](https://conda.io/docs/user-guide/install/index.html
Assuming that you have a 64-bit system, on Linux, download and install Anaconda 3
```
$ wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
$ bash Anaconda3-5.0.1-Linux-x86_64.sh
```
Answer yes to the user agreement; go with the default installation path or specify your own. Answer yes to prepend the install location to your PATH.

# Step 2: Prepare analysis directory
## Create a new directory and change into it.
## Download data files from the tar archive - this will set up a data/ folder with annotation files used in for the manuscript
## Clone this repository

# Step 3: Create environment `regulatory_comparisons` with the required software using the `environment.yaml` file and activate it
```
$ conda env create --name regulatory_comparisons --file environment.yaml
$ source activate regulatory_comparisons
```
Now you can use the installed tools. Check if all went fine:
```
$ snakemake --help
```
	 
