import pandas

pandas.options.mode.chained_assignment = None

print(snakemake.input.geneList[0])
""" Read ESI results computed previously """
dGeneList = pandas.read_csv(snakemake.input.geneList[0], sep='\t', usecols=['Name','Description',snakemake.params.cellName])
dGeneList.columns = ['Name','Description','lclesi']

""" Read Gencode V19 TSS file """
dgencode_annotation = pandas.read_csv(snakemake.input.gencode_annotation, sep='\t', header=None, names=["chrom","start","end","Name","Description","type"], usecols=["chrom","start","end","Name","Description"])

d = pandas.merge(dgencode_annotation, dGeneList, how="inner", on=['Name','Description'])
d[['chrom','start','end','Name']].to_csv(snakemake.output.fileFullGeneSet, sep='\t', index=False, header=False)



