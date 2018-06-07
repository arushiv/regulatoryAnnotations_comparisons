# Use `Snakefile` to compute:
- Number of segments
- Segment length distribution
- Percent genome coverage
and plot these summary statistics for annotations (Fig. 1 A,B,C)
- compute overlap fractions and plot heatmap (Fig. 1D)
Does not require cluster submission.

## Dry run
```
$	snakemake -np
```
## Run
```
$	snakemake -p --latency-wait 400
```
	
## Print workflow DAG
```
	snakemake --dag | dot -Tsvg > workflow.svg
```

# Use `Snakefile_GAT` to compute:
- Overlap enrichment between each pair of regulatory annotations. (Supplementary Fig 1)

## Dry run
```
$	snakemake -nps Snakefile_GAT
```
## Run
```
$ mkdir -p logs && snakemake --cluster-config cluster.yaml  --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --job-name {cluster.jobname} -o {cluster.output} -e {cluster.error} --parsable" -j 60 --latency-wait 400 -s Snakefile_GAT --configfile config_GAT.yaml  -p
```
	
## Print workflow DAG
```
	snakemake -s Snakefile_GAT --dag | dot -Tsvg > workflow_GAT.svg
```
