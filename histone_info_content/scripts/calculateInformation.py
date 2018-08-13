import pandas
import numpy
import os
import sys
import math
import argparse

def log_with_nan(x, y):
    try:
        return math.log(x, y)
    except ValueError:
        """Places with log(0)"""
        return float('nan')

def entropy_function(x):
    y = -x * log_with_nan(x, 2)
    return y

def calculate_information(d):
    """
    Relative RPKM = x(g,t)/sum(x(g,t))
    Entropy = - sum(relativeRPKM(g,t) * log2(relative RPKM(g,t)))
    I = Relative RPKM * (Entropy max - Entropy(g) )
    Entropy max (all probabilities equal)= sum(-1/N * log2(-/N)) = log2(N)
    """

    N = len(d.columns)
    print(f"Number of features = {N}")
    
    entropy_max = math.log(N, 2)
    print(f"max entropy = {entropy_max}")
    
    d_relative_rpkm = d.div(d.sum(axis=1), axis=0)
    d_relative_rpkm.loc[:,'entropy'] = d_relative_rpkm.applymap(lambda x: entropy_function(x)).sum(axis=1)

    d_information = d_relative_rpkm.apply(lambda x: x*(entropy_max - x['entropy']), axis=1)
    d_information.drop('entropy', axis=1, inplace=True)
    d_information.loc[:,'change_in_entropy'] = entropy_max - d_relative_rpkm['entropy']
    return d_information




def getOpts():
    parser = argparse.ArgumentParser(description='For a given matrix of say, RPKM values for different GTEx samples, calculate mean RPKM for each tissue.'
                                     'Then for each row, calculate information content and information for each column',
                                     usage='python calculateInformation.py <inputfile> <outputfile>')
    parser.add_argument('datafile', help="""Input Matrix. Tab separated, with header specifying sample IDs or tissues.""")
    parser.add_argument('outputfile', help="""Output file.""")
    parser.add_argument('--indexCols', nargs='+', help="""Columns other than those across which information content has to be calculated.""")
    parser.add_argument('--subsetCols', nargs='+', help="""Columns to subset on before calculating information. Rest columns will be left out.""")
    parser.add_argument('--cell', type=str, help="""Columns name which would be copied to the 'information' column to concat and plot later.""")
    parser.add_argument('-e', '--expressionFilter', type=float, help = """Min expression for value to filter for. """)
    parser.add_argument('-et', '--expressionFilterTissue', type=str, help = """Tissue name for which genes have to be filtered based on gene expression in the 'expressionFilter' flag""")

    return parser


def runprog():
    t = pandas.read_csv(args.datafile, sep='\t')
    t.loc[:,'avg_post'] = t['avg_posterior']
    t.set_index(['chrom', 'start', 'end', 'infoContent', 'cell', 'annotation', 'avg_post'], inplace=True)
    t.rename(columns={'avg_posterior':args.cell}, inplace=True)

    if args.subsetCols is not None:
        print(args.subsetCols)
        t = t[args.subsetCols]


    """
    Calculate information
    """
    out = calculate_information(t)
    out.reset_index(inplace=True)
    out.loc[:,'information'] = out[args.cell]
    return out


if __name__ == '__main__':
    parser = getOpts()
    args = parser.parse_args()

    # inputdf = pandas.read_csv(args.datafile, sep='\t', index_col=args.indexCols)
        
    # if (args.expressionFilter is not None) and (args.expressionFilterTissue is not None):
    #     inputdf = filterGeneExpression(inputdf, args.expressionFilter, args.expressionFilterTissue)
    #     print(len(inputdf.index))

    try:
        out = runprog()
    except pandas.errors.EmptyDataError:
        out = pandas.DataFrame()    
    
    out.to_csv(args.outputfile, sep='\t', na_rep="NA", index=False)
