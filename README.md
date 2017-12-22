# Regulatory Annotations Comparisons
This repository contains code to reproduce results of the manuscript:

Each directory contains code for analyses specific to different figures of the paper. Everything is run using `snakemake`. These directories follow this general pattern:
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

# Preliminaries

- Install [Anaconda](https://conda.io/docs/user-guide/install/index.html)
- Install [snakemake](http://snakemake.readthedocs.io/en/latest/getting_started/installation.html) (This code used version 4.2.0)


	 
