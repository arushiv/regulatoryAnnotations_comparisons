Compute -
- Number of segments
- Segment length distribution
- Percent genome coverage

and plot these summary statistics for annotations (Fig. 1 A,B,C)

Also compute overlap fractions and plot heatmap (Fig. 1D)
Does not require cluster submission.

.. image:: https://github.com/arushiv/regulatoryAnnotations_comparisons/blob/master/summaryStatistics/workflow.svg

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