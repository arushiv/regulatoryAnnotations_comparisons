Analyses for Fig. 5A, B 
Plot effect size distribution of GTEx V7 Whole Blood eQTL in annotations and run power calculation. Also run regression model for supplementary table 1
	
## Dry run
```
$	snakemake -np --configfile config_files/config.yaml
```
	
## Print rule graph
```
$ snakemake --configfile config_files/config.yaml --rulegraph | dot -Tsvg > workflow.svg
```

## Run for SLURM: 
```
$	mkdir -p logs && snakemake --cluster-config cluster.yaml  --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --job-name {cluster.jobname} -o {cluster.output} -e {cluster.error} --parsable" -j 60 -p --latency-wait 400  --configfile config_files/config.yaml
```
	