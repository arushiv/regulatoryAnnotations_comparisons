import pandas

desi = pandas.read_csv(snakemake.input.esi, sep='\t', usecols=['Name', 'Description', snakemake.params.cellName])
desi.columns = ['Name','Description','lclesi']

deqtl = pandas.read_csv(snakemake.input.eqtl, sep='\t', usecols=['snp', 'GENE_ID', 'rvalue', 'Position', 'Alleles'])
deqtl.loc[:,'Name'] = deqtl['GENE_ID'].map(lambda x: x.split('.')[0])
deqtl.drop(['GENE_ID'], inplace=True, axis=1)

d = pandas.merge(desi, deqtl, how="inner", on=['Name'])
d.loc[:,'lclesi_bin'] = pandas.qcut(d['lclesi'], q=snakemake.params.bins, labels=list(range(1, snakemake.params.bins + 1)))

d.to_csv(snakemake.output.fullfile, sep='\t', index=False)
