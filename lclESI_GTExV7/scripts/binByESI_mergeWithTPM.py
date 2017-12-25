import pandas
import numpy

esi = pandas.read_csv(snakemake.input.esi[0], sep='\t', usecols = ['Name','Description','Cells - EBV-transformed lymphocytes'])
esi.columns = ['Name','Description','lclesi']

""" Bin by LCL ESI """
esi.loc[:,'lclesi_bin'] = pandas.qcut(esi['lclesi'], q=5, labels=list(range(1, 6)))

""" Merge with RPKM info """
rpkm = pandas.read_csv(snakemake.input.rpkm, sep='\t', comment="#", header=1)
rpkm.rename(columns={'gene_id':'Name'}, inplace=True)
rpkm.loc[:,'Name'] = rpkm['Name'].map(lambda x: x.split('.')[0])
# rpkm = rpkm[rpkm['Cells - EBV-transformed lymphocytes'] >= float(snakemake.params.value)]
d = pandas.merge(esi, rpkm, how="inner", on=['Name','Description'])
# d = d[d['Cells - EBV-transformed lymphocytes'] >= float(snakemake.params.value)]
# d.loc[:,'lclesi_bin'] = pandas.qcut(d['lclesi'], q=5, labels=False)

drpkm = pandas.melt(d, id_vars=['Name','Description', 'lclesi', 'lclesi_bin'], var_name='cell', value_name='tpm')
drpkm.to_csv(snakemake.output[0], sep='\t', index=False, na_rep="NA")  
