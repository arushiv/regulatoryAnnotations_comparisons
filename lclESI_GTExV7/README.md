## Run
```
$       snakemake -p --latency-wait 400
```
	
Use GTEx v7 RNA seq data to calculate LCL-ESI
This workflow is included as a sub-workflow of analyses that require these results, so no need to run this separately.
	
Note:
For genes with very low expression across tissues, low ESI may be over-represented. Filter for minimum expression value and check if the lclESI distribution looks roughly normal.
After testing a number of minimum LCL TPM filters, min LCL median TPM = 0.15 was selected.	
The ESI values for these set of genes were used in all in further analyses.

