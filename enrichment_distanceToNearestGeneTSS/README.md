# Use `Snakefile_allPCgenes` to:
Calculate cumulative distribution of distance to nearest protein-coding gene TSS for each regulatory annotation (Supplementary Fig 3)

## Run:
```
snakemake -ps Snakefile_allPCgenes --latency-wait 400
```

# Use `Snakefile_binBylclESI` to:
Split genes into bins my lclESI (calculated in a subworkflow). Get TSS for genes in each bin. Calculate enrichment as the cumulative distribution of distance to nearest gene in each bin, divided by cumulative distribution of distance to nearest gene from a randomly sampled, equal-sized gene set from across bins (Fig. 4A)

## Run for enrichment Fig. 4A: cluster execution recommended. Example for SLURM: 
```
$	mkdir -p logs && snakemake --cluster-config cluster.yaml  --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --job-name {cluster.jobname} -o {cluster.output} -e {cluster.error} --parsable" -j 60 -p --latency-wait 400 -s Snakefile_binBylclESI
```
	