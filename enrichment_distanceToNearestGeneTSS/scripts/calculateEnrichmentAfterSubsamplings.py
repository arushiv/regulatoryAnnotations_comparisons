#!/usr/bin/env python

import argparse
import glob
import os
import pandas
import pybedtools
import numpy
from statsmodels.distributions.empirical_distribution import ECDF
import scipy

"""
Script function:
Compute enrichment for distance to nearest gene
"""

def ecdfOfDataframeColumn(x):
    ecdf = ECDF(x)
    rangeOfx = numpy.arange(0, 8.5, step=0.05)
    return ecdf(rangeOfx)

def getname(x):
    return os.path.basename(x)

def getXColumn(x):
    x.loc[:,'xvals'] = numpy.arange(0, 8.5, step=0.05).tolist()
    return x
    
def getEcdf(df, groupings):
    df = df[df['distance'] != -1]   # pybedtools.closest returns -1 if no feature is found on the same chromosome. Remove this from the output dataframe
    df.loc[:,'distance'] = numpy.log10(df['distance'] + 1)
    fulldf = df.groupby(groupings)['distance'].apply(pandas.core.groupby.Series.tolist).apply(ecdfOfDataframeColumn).reset_index()
    fulldf = fulldf.rename(columns={'distance': 'ecdf_y'})
    fulldf = fulldf['ecdf_y'].apply(pandas.Series).stack().reset_index(level=1, drop=True).to_frame('ecdf_y').join(fulldf[groupings], how='left')
    fulldf = fulldf.groupby(groupings).apply(getXColumn)
    return fulldf


def getMean(df, groupings):
    outdf = df.groupby(groupings).agg([numpy.mean, scipy.stats.sem]).reset_index()
    # outdf = df.groupby(groupings).agg([numpy.mean, numpy.std(ddof=1)]).reset_index()   # standard deviation

    outdf.columns = outdf.columns.droplevel()
    groupings.extend(['ecdf_y_shuffleMean','ecdf_y_shuffleSem'])
    outdf.columns = groupings
    return outdf

def printOutput(outdf, outfile):
    outdf.to_csv(outfile, index=False, sep='\t', header=True, na_rep="NA")

def getOpts():
    parser = argparse.ArgumentParser(description="""Compute enrichment for nearest distance to TSS with genes.""")
    parser.add_argument('mainDistanceAnnotFile', help="""Data file with nearest gene TSS distance values for each regulatory annotation for which enrichment has to be calculated.'""")
    parser.add_argument('outfile', type=str, default="outfileshuffle.dat", help="""Output file name.""")
    parser.add_argument('-e', '--ecdfFileList', nargs='+', type=str, help="""List of multithreaded ECDF files of random subsamplings""")
    return parser
    
if __name__ == '__main__':

    parser = getOpts()
    args = parser.parse_args()
    pandas.options.mode.chained_assignment = None
   
    """
    Get ECDF for dex response gene nearest distance
    """
    distanceAnnotFile = pandas.read_csv(args.mainDistanceAnnotFile, sep='\t')
    groupings_distAnnot = ['cell','annotation','bin']
    distanceAnnotFile_ecdf = getEcdf(distanceAnnotFile, groupings_distAnnot)   #.to_csv("test.df", sep='\t', index=False)
    distanceAnnotFile_ecdf = distanceAnnotFile_ecdf.round({'xvals' : 2})
    
    """
    Concat Shuffling ECDFs
    """
    def readECDFs(ecdfFile):
        return pandas.read_csv(ecdfFile, sep='\t')
    
    shuffledf = pandas.concat([readECDFs(ecdfFile) for ecdfFile in args.ecdfFileList])
    shuffledf = shuffledf.round({'xvals' : 2})

    """    
    Compute mean and std dev of ECDF from shufflings
    """
    groupings_shuffleMeans = ['cell','annotation','xvals']
    shuffle_MeanEcdf = getMean(shuffledf, groupings_shuffleMeans)

    outdf = pandas.merge(distanceAnnotFile_ecdf, shuffle_MeanEcdf, on=['cell','annotation','xvals'])

    printOutput(outdf, args.outfile)
