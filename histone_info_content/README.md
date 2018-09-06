# Use `Snakefile` to:
- Chromatin state enhancer and promoter posterior probability and information content analysis for Fig 3 A and B.

## Dry run
```
$	snakemake -np
```
## Run
```
$	snakemake -p --latency-wait 400
```
	
## Print workflow rulegraph
```
	snakemake --rulegraph | dot -Tsvg > workflow.svg
```

## Run
```
$ mkdir -p logs && snakemake --cluster-config cluster.yaml  --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --job-name {cluster.jobname} -o {cluster.output} -e {cluster.error} --parsable" -j 60 --latency-wait 400  --configfile  config_files/config_chromhmm.yaml -p
```
	
