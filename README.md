# Regulatory Annotations Comparisons
This repository contains code to reproduce results of the manuscript:

Each directory contains code and dependency specifications for analyses specific to different figures of the paper. These directories follow this general pattern:
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

- Install [snakemake](https://bitbucket.org/johanneskoester/snakemake) (This code used version 4.2.0)


	 
