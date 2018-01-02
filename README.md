# Regulatory Annotations Comparisons
This repository contains code to reproduce results of the manuscript:

Each directory contains code for analyses specific to different figures of the paper. Everything is run using [snakemake](http://snakemake.readthedocs.io/en/latest/). The analysis directories follow this general pattern:
```	
├── .gitignore
├── README.md
├── config.yml : Cluster config specifications (if analysis is recommended to be run on a cluster)
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
* bedtools 2.26.0
* GREGOR 1.2.1
* GAT 1.3.5

To setup these pre-requisites, use the Anaconda/Miniconda Python3 distribution. The Conda package manager is used to obtain and deploy the defined software packages in the specified versions. These instructions are for the Linux platform
	
# Step 1: Install [Anaconda3](https://conda.io/docs/user-guide/install/index.html)
Assuming that you have a 64-bit system, on Linux, download and install Anaconda 3
```
$ wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
$ bash Anaconda3-5.0.1-Linux-x86_64.sh
```
Answer yes to the user agreement; go with the default installation path or specify your own. Answer yes to prepend the install location to your PATH.

Required for analyses computing enrichment of eQTL and GWAS in regulatory annotations, please install [GREGOR](https://genome.sph.umich.edu/wiki/GREGOR) manually.

# Step 2: Prepare analysis directory
Create a new directory and change into it.
Clone this repository. This will produce the directory structure to replicate specific analyses. Edit the `BASE_PATH` in the `Snakefile_config`

# Step 3: Create and activate environment
Create an environment named `regulatory_comparisons` with the required software using the `environment.yaml` file and activate it
```
$ conda env create --name regulatory_comparisons --file environment.yaml
$ source activate regulatory_comparisons
```
Now you can use the installed tools. Check if all went fine:
```
$ snakemake --help
```
# Step 4: Execute analyses
Now change into the analysis directories and run snakemake using the respective Snakefiles. Each Snakefile will check for and download required data if missing as per the `Snakefile_config` in the top directory.

Below are the analysis directories corresponding to each figure:
```	
├── Fig 2A-D: `summaryStatistics` :: `Snakefile`
├── Fig 3: `coverageChromatinStates` :: `Snakefile`
├── Fig 4A: `enrichment_distanceToNearestGeneTSS` :: `Snakefile_binBylclESI`
│   ├── Supplementary Fig. S2 :: `Snakefile_allPCgenes`	
│   ├── Supplementary Fig S4: Automatically runs subworkflow in `lclESI_GTExV7`
│   └──	Supplementary Fig S5 :: `Snakefile_binBylclESI`
├── Fig 4B: `enrichment_eQTL` :: `Snakefile_binByESI`
│   ├── Supplementary Fig. S3 :: `Snakefile_bulkEnrichment`	
│   ├── Automatically runs subworkflow in `lclESI_GTExV7`
│   └──	Supplementary Fig S6 :: `Snakefile_binByESI`	
├── Fig 5: `qtl_effectSizeDistribution`
    ├── Fig 5A :: `Snakefile_eQTLEffect`
    ├── Supplementary Fig S7 :: `Snakefile_eQTLEffect`	
    ├── Fig 5B :: `Snakefile_dsQTLEffect`
    └── Fig 5C :: `Snakefile_allelicBiasEffect`

```
	