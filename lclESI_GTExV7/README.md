Use GTEx v7 RNA seq data to calculate LCL-ESI
GTEx v7 provides expression quantifications in TPM.
Note:
For genes with very low expression across tissues, low ESI may be over-represented. Filter for minimum expression value and check if the lclESI distribution looks roughly normal in `Snakefile_testTPMThresholds`
After testing a number of minimum LCL TPM filters, min LCL median TPM = 0.15 was selected.	
The ESI values for these set of genes were used in all in further analyses.