#!/usr/bin/env python
import argparse
import pandas
import math


def calculateMean(datafile, samplefile):
        sampledf = samplefile.groupby('tissue')['sample'].apply(list)

        """Take those tissues that have more than 25 samples"""

        sampledf = sampledf[sampledf.apply(lambda x: len(x)) > 25]   
        inputdf = pandas.concat(sampledf.apply(lambda x: chunk.loc[:,x].mean(axis=1)).transpose() for chunk in datafile)
        
        return inputdf

def log_with_nan(x, y):
        try:
                return math.log(x, y)
        except ValueError:
                """Places with log(0)"""
                return float('nan')

def entropy_function(x):
        y = -x * log_with_nan(x, 2)
        return y

def filterProteinCoding(d, geneFilterList, filterName):
        dCoding = pandas.read_csv(geneFilterList, sep='\t')
        dCoding = dCoding[dCoding['type'] == filterName]

        """
        Rename Gene ID to remove '.' to match with gene list
        """
        d.reset_index(inplace=True)
        """ Column name is gene_id in GTEx v7 - convert to 'Name' to be consistent with GTEx v6p """
        d.rename(columns={'gene_id' : 'Name'}, inplace=True)
        d.loc[:,'Name'] = d['Name'].map(lambda x: x.split(".")[0])
        outdf = pandas.merge(d, dCoding[['Name']], how="inner", on=['Name'])
        outdf.set_index(['Name','Description'], inplace=True)
        return outdf

def filterGeneExpression(d, expressionFilter, expressionFilterTissue):
        return d[d[expressionFilterTissue] >= expressionFilter]

def calculate_ESI(d):
        """
        Relative RPKM = x(g,t)/sum(x(g,t))
        Entropy = - sum(relativeRPKM(g,t) * log2(relative RPKM(g,t)))
        Q(g,t) = Entropy(g) - log2(relativeRPKM(g,t))
        maxQ = max(Q(t))
        """
        d_relative_rpkm = d.div(d.sum(axis=1), axis=0)  
        d.loc[:,'entropy'] = d_relative_rpkm.applymap(lambda x: entropy_function(x)).sum(axis=1)   
        d_qvalues = d_relative_rpkm.applymap(lambda x: log_with_nan(x, 2)).apply(lambda x: d.loc[:,'entropy'] - x, axis=0)
        # print(d_qvalues.loc[d_qvalues[['Cells - EBV-transformed lymphocytes']].idxmax(), 'Cells - EBV-transformed lymphocytes'])
        d_max_q = d_qvalues.max(axis=0).to_frame().transpose()       
        d_esi = d_qvalues.apply(lambda x: 1 - x/d_max_q.squeeze(), axis=1)
        return d_esi

def getOpts():
        parser = argparse.ArgumentParser(description='For a given matrix of say, RPKM values for different GTEx samples, calculate mean RPKM for each tissue. Sample IDs corresponding each tissue are specified in another file. Option to skip calculating mean if sample ID file is not provided, matrix file headers are then taken to be tissue types. Then calculate ESI for each tissue', usage='python esiScoreAfterMean.py datafile.txt -s samplefile.txt output.txt')
        parser.add_argument('datafile', help="""Input Matrix. Tab separated, with header specifying sample IDs or tissues. First two columns taken as 'Name' and 'Description'.""")
        parser.add_argument('-s', '--samplefile', help="""If mean for each tissue has to be calculated, supply a tab separated file with 2 columns with headers: 'tissue' 'sample', specifying sample IDs corresponding to each tissue type.""")
        parser.add_argument('outputfile', help ="""Output file.""")
        parser.add_argument('-f', '--geneFilterList', help = """File with gene 'Name' as first column and gene 'type' as second column. No header. will filter for genes with type = "protein_coding" before calculating ESI""")
        parser.add_argument('-t', '--geneFilterType', help = """value of the 'type' column on the geneFilterList on which genes should be filter. Example, give 'protein_coding' to select these genes from the provided list before calculating ESI""")
        parser.add_argument('-e', '--expressionFilter', type=float, help = """Min expression for value to filter for. """)
        parser.add_argument('-et', '--expressionFilterTissue', type=str, help = """Tissue name for which genes have to be filtered based on gene expression in the 'expressionFilter' flag""")

        return parser
        
if __name__ == '__main__':
        parser = getOpts()
        args = parser.parse_args()
        
        outputfile = args.outputfile
        
        """
        If a sample vs tissue file is provided, calculate Mean RPKM for each tissue and then proceed
        """
        
        if args.samplefile is not None:
                """
                Read input matrix - expected to be very large, so is read in chunks. Note first line is skipped while reading header because of the specific GTEx input file formatting.
                Firt two columns - gene name and description are read as index. 
                """
                datafile = pandas.read_csv(args.datafile, sep='\t', index_col=[0,1], comment='#', header=1, low_memory=False, chunksize = 10000, iterator=True)
                samplefile = pandas.read_csv(args.samplefile, sep='\t')
                inputdf = calculateMean(datafile, samplefile)

        else:
                """
                If using the GTEx V6p median gene RPKM per tissue file, no need to calculate mean
                """
                inputdf = pandas.read_csv(args.datafile, sep='\t', index_col=[0,1], comment='#', header=1)

        if args.geneFilterList is not None:
                inputdf = filterProteinCoding(inputdf, args.geneFilterList, args.geneFilterType)
        print(len(inputdf.index))

        if (args.expressionFilter is not None) and (args.expressionFilterTissue is not None):
                inputdf = filterGeneExpression(inputdf, args.expressionFilter, args.expressionFilterTissue)
        print(len(inputdf.index))

        """
        Calculate ESI
        """
        d_esi = calculate_ESI(inputdf)
        d_esi.to_csv(outputfile, sep='\t', na_rep="NA")
