#!/usr/bin/env python
import argparse
import glob
import os
import pandas
import pybedtools
import numpy
from statsmodels.distributions.empirical_distribution import ECDF
import subprocess as sp

def fixFunc(df, bedfile):
    # ndf = df.closest(pybedtools.BedTool(bedfile), d=True).to_dataframe()  # Closest for genes
    ndf = pybedtools.BedTool(bedfile).closest(df, d=True).to_dataframe()  # Closest for annotations
    ndf = ndf[ndf.columns[-1:]]  # Retain only last column of the closest distance
    ndf.columns = ['distance']
    ndf.loc[:,'name'] = os.path.basename(bedfile).replace(".annotations.bed", "")
    return ndf

def getNearestDistance(df, bedFileList):
    returndf = pandas.concat([fixFunc(df, bedfile) for bedfile in bedFileList])
    returndf[['cell','annotation']] = returndf['name'].str.split('.', expand=True)
    returndf.drop(['name'], axis=1)
    return returndf

def ecdfOfDataframeColumn(x):
    ecdf = ECDF(x)
    rangeOfx = numpy.arange(0, 8.5, step=0.05)
    return ecdf(rangeOfx)

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


def printOutput(df, outfile):
    df.to_csv(outfile, sep='\t', index=False, na_rep="NA")

def getDistDataframe(distanceFilesList):
    listOfDFs = [pandas.read_csv(distanceFile, sep='\t', header=None, names=['distance','cell','annotation']) for distanceFile in distanceFilesList]
    distDf = pandas.concat(listOfDFs)
    return distDf

def getOpts():
    parser = argparse.ArgumentParser(description="""Shuffling: Given a  bed file, randomly sample on each unique item in name or other column. 2. Go back to original bedfile and fetch all coordinates corresponding to sampled name field. 3. Compute bedtools closest given a list of annotation files for each sampling iteration and return dataframe.""")
    parser.add_argument('--tssFile', type=str, help="""bed file with gene TSS coordinates. No header""")
    parser.add_argument('--bedFileList', nargs='+', help="""Annotation bed files to get nearest dist and ECDF.""")
    parser.add_argument('--distanceFiles', nargs='+', help="""File with ditance to nearest gene TSS if already calculated. No need for --bedFileList and --tssFile in this case""")
    parser.add_argument('--outputfile', required=True, type=str, help="""Output file name.""")

    return parser

    
if __name__ == '__main__':

    parser = getOpts()
    args = parser.parse_args()

    pandas.options.mode.chained_assignment = None
    
    if args.bedFileList is not None:
        """
        Get Nearest Distance for annotations
        """
        
        tssFile = pybedtools.BedTool(args.tssFile).sort()
        distDf = getNearestDistance(tssFile, args.bedFileList)

    elif args.distanceFiles is not None:

        distDf = getDistDataframe(args.distanceFiles)

    else:
        print("Provide list of annotation bed files or Distance to nearest gene files")
    
    """    
    Get ECDF for shufflings
    """

    groupings = ['cell','annotation']
    df_Ecdf = getEcdf(distDf, groupings)
    printOutput(df_Ecdf, args.outputfile)
        
    


