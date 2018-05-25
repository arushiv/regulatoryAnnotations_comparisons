Analyses for Fig. 5C 
Plot effect size distribution of DNase QTL in annotations
## Dry run
```
$	snakemake -np
```
## Run for SLURM: 
```
$	mkdir -p logs && snakemake --cluster-config cluster.yaml  --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --job-name {cluster.jobname} -o {cluster.output} -e {cluster.error} --parsable" -j 60 -p --latency-wait 400 
```
	