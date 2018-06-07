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
    ndf = pybedtools.BedTool(bedfile).closest(df, d=True, t="first").to_dataframe()  # Closest for annotations, report one for ties
    ndf = ndf[ndf.columns[-1:]]  # Retain only last column of the closest distance
    ndf.loc[:,'name'] = bedfile
    return ndf

def sampleUnique(tssFile, colNumber, numberOfGenes):
    genelist = tssFile.iloc[:,colNumber].drop_duplicates().sample(numberOfGenes)
    return pybedtools.BedTool.from_dataframe(tssFile[tssFile.iloc[:,colNumber].isin(genelist)]).sort()

def getNearestDistance(df, bedFileList, runNumber):
    returndf = pandas.concat([fixFunc(df, bedfile) for bedfile in bedFileList])
    returndf.loc[:,'runNumber'] = runNumber
    return returndf

def getname(x):
    return os.path.basename(x)

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

def getOpts():
    parser = argparse.ArgumentParser(description="""Shuffling: Given a  bed file, randomly sample on each unique item in name or other column. 2. Go back to original bedfile and fetch all coordinates corresponding to sampled name field. 3. Compute bedtools closest given a list of annotation files for each sampling iteration and return dataframe.""")
    parser.add_argument('tssFile', type=str, help="""bed file on which sampling has to be done. No header""")
    parser.add_argument('outfile', type=str, default="outfileshuffle.dat", help="""Output file name.""")

    parser.add_argument('-b', '--bedFileList', nargs='+', type=str, help="""List of annotation bedfiles [3 column tab separated no header] to compute closest distance after each sampling""")
    parser.add_argument('-s', '--numberOfGenesFile', type=str, help="""File with gene names that were considered to calculate ESI""")
    parser.add_argument('--bins', type=int, help="""Number of ESI Bins """)
    parser.add_argument('-c', '--colNumber', type=int, default=4, help="""Column number of gene names from which unique genes have to be sampled""")
    parser.add_argument('-r', '--runs', type=int, default=5, help="""Number of sampling runs to perform. Data will be recorded with additional column of the run number. Default = 5 runs""")
    parser.add_argument('-t', '--threadNumber', type=int, help="""If threading 10 jobs in Snakemake, specify which thread to save in dataframe. This will be appended to the run number column so each subsampling of each thread can be accounted for""")
    return parser

def getNumberOfGenes(filename, bins):
    cmd = r"""less {filename} | grep -v Name | cut -f1 | sort | uniq | wc -l """.format(filename=filename)
    n=int(sp.check_output(cmd, shell=True).decode("utf-8").rstrip())
    return(int(n/bins))
    
if __name__ == '__main__':

    parser = getOpts()
    args = parser.parse_args()

    pandas.options.mode.chained_assignment = None
    
    bedFileList = args.bedFileList

    if bedFileList is not None:
        tssFile = pandas.read_csv(args.tssFile, header=None, sep='\t')
        numberOfGenes = getNumberOfGenes(args.numberOfGenesFile, args.bins)
        print(numberOfGenes)
        colNumber = args.colNumber - 1   # python index from 0
        
        """
        Shufflings
        """
        def samplingResult(runNumber):
            print(runNumber)
            tdf = sampleUnique(tssFile, colNumber, numberOfGenes)
            distDf = getNearestDistance(tdf, bedFileList, runNumber)
            return distDf

        RUNS = ['s{0}_{1}'.format(run, args.threadNumber) for run in list(range(args.runs))]    # Make string list of all shuffle runs to perform
        shuffledf = pandas.concat([samplingResult(runNumber) for runNumber in RUNS])
        shuffledf.columns = ['distance', 'name', 'runNumber']

        shuffledf.loc[:,'cell'], shuffledf.loc[:,'annotation'], shuffledf.loc[:,'extra'] = shuffledf.loc[:,'name'].str.split('.', 2).str
        shuffledf.loc[:,'cell'] = shuffledf['cell'].map(os.path.basename)
        shuffledf = shuffledf[['distance','cell','annotation','runNumber']]
    
        """    
        Get ECDF for shufflings
        """
        groupings_shuffle = ['cell','annotation','runNumber']
        shuffle_Ecdf = getEcdf(shuffledf, groupings_shuffle)
        
        printOutput(shuffle_Ecdf, args.outfile)
        
    else:
        print("Provide list of annotation bed files")
    


