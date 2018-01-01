- `Snakefile_bulkEnrichment` : Calculate enrichment of GTEx LCL in annotations (Supplementary fig S3)
- `Snakefile_binByESI` : Bin eQTL by lclESI and calculate enrichment of eQTL sets. (Fig. 4B)

## Dry run
```
$	snakemake -nps Snakefile_binByESI
```
## Run for supplementary figure S3: cluster execution recommended. Example for SLURM:
```
$	snakemake --cluster-config cluster.yml  --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --job-name {cluster.jobname} -o {cluster.output} -e {cluster.error} --parsable" -j 60 -p --latency-wait 400 -s Snakefile_bulkEnrichment
```
## Run for enrichment Fig. 4B: cluster execution recommended. Example for SLURM:
```
$	snakemake --cluster-config cluster.yml  --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --job-name {cluster.jobname} -o {cluster.output} -e {cluster.error} --parsable" -j 60 -p --latency-wait 400 -s Snakefile_binByESI
```
	