# dry run
snakemake -np

# Run
snakemake -p --latency-wait 400

# print workflow
snakemake --dag | dot -Tsvg > workflow.svg
