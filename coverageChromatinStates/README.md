# Use `Snakefile` to compute:
	- Coverage of annotations in chromatin states
	A subworkflow determines and downloads data files required for this analysis.
	For each annotation, plot fraction of overlap in chromatin states across 4 cell types

Cluster execution is not required
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
$	snakemake --dag | dot -Tsvg > workflow.svg
```
	