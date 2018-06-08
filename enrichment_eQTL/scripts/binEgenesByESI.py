import pandas

desi = pandas.read_csv(snakemake.input.esi, sep='\t', usecols=['Name', 'Description', snakemake.params.cellName])
desi.columns = ['Name','Description','lclesi']

deqtl = pandas.read_csv(snakemake.input.eqtl, sep='\t', usecols=['gene_id', 'gene_name', 'slope', 'chrom', 'pos', 'ref', 'alt'])
deqtl.loc[:,'Name'] = deqtl['gene_id'].map(lambda x: x.split('.')[0])
deqtl.drop(['gene_id'], inplace=True, axis=1)

d = pandas.merge(desi, deqtl, how="inner", on=['Name'])
d.loc[:,'lclesi_bin'] = pandas.qcut(d['lclesi'], q=snakemake.params.bins, labels=list(range(1, snakemake.params.bins + 1)))

for name, group in d.groupby('lclesi_bin'):
    filename = "{namestring}{name}.txt".format(namestring=snakemake.params.namestring, name=name)
    group[['chrom','pos']].drop_duplicates().to_csv(filename, index=False, header=False, sep=':')

d.to_csv(snakemake.output.fullfile, sep='\t', index=False)
