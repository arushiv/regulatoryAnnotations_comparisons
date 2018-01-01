Analyses for Fig. 4A and Supplementary Fig. S2 
Split genes into bins my lclESI (calculated in a subworkflow)
Get TSS for genes in each bin. Calculate enrichment as the cumulative distribution of distance to nearest gene in each bin, divided by cumulative distribution of distance to nearest gene from a randomly sampled, equal-sized gene set from across bins

## Dry run
```
$	snakemake -nps Snakefile_binBylclESI
```
## Run for supplementary figure
```
snakemake -ps Snakefile_allPCgenes --latency-wait 400
```
## Run for enrichment Fig. 4A: cluster execution recommended. Example for SLURM: 
```
$	snakemake --cluster-config cluster.yml  --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --job-name {cluster.jobname} -o {cluster.output} -e {cluster.error} --parsable" -j 60 -p --latency-wait 400 -s Snakefile_binBylclESI
```
	