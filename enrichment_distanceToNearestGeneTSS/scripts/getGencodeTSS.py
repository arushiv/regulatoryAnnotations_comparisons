import pandas

pandas.options.mode.chained_assignment = None

""" Read ESI results computed previously """
dGeneList = pandas.read_csv(snakemake.input.geneList, sep='\t', usecols=['Name','Description',snakemake.params.cellName])
dGeneList.columns = ['Name','Description','lclesi']

""" bin by ESI """
dGeneList.loc[:,'lclesi_bin'] = pandas.qcut(dGeneList['lclesi'], q=snakemake.params.bins, labels=list(range(1, (snakemake.params.bins + 1))))

""" Read Gencode V19 TSS file """
dgencode_annotation = pandas.read_csv(snakemake.input.gencode_annotation, sep='\t', header=None, names=["chrom","start","end","Name","Description","type"], usecols=["chrom","start","end","Name","Description"])

d = pandas.merge(dgencode_annotation, dGeneList, how="inner", on=['Name','Description'])
d[['chrom','start','end','Name']].to_csv(snakemake.output.fileFullGeneSet, sep='\t', index=False, header=False)

""" Make separate TSS files for each bin """
for name, group in dGeneList.groupby(['lclesi_bin']):
    df = pandas.merge(dgencode_annotation, group, how="inner", on=['Name','Description'])
    filename = "{filenameString}{name}.bed".format(filenameString=snakemake.params.filenameString, name=name)
    df.to_csv(filename, sep='\t', index=False)


