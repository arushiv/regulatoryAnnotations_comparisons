# Run GREGOR: Fix conf files, submit drmr commands, compile dataframe and plot 
#!/usr/bin/env python

import os
import argparse
import pandas
import glob


def makedataframe(filename, fieldSeparator, outputfilename):
    framelist = []
    for name in filename:
        job_info = os.path.basename(name.replace('/StatisticSummaryFile.txt','')).replace('output_','')

        df = pandas.read_csv(name, sep='\t')
        df['Bed_File'] = df['Bed_File'].str.replace('.bed','')
        df['job_info'] = job_info
        splitBeddf = df['Bed_File'].str.split(fieldSeparator, expand=True)
        bincol = df['job_info'].str.split("_", expand=True)[1].str.replace('lclESIbin', '')
        # splitSnpdf = splitJobdf.iloc[:,0].str.split(".", expand=True)
        # outdf = pandas.concat([splitSnpdf, splitJobdf.drop(splitJobdf.columns[0], 1), splitBeddf, df.drop(['Bed_File','job_info'], 1)], 1)
        odf = pandas.concat([bincol, splitBeddf, df.drop(['Bed_File'], 1)], 1)
        framelist.append(odf)
    outdf = pandas.concat(framelist)
    return outdf
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compile dataframe from GREGOR output folders.')
    parser.add_argument('--filename', nargs='+', type=str, help="""The StatisticSummaryFile from GREGOR.""")
    parser.add_argument('-f','--nameFieldSeparator', type=str, default='.', help="""Field separator to make columns from bed file name. (Default='.')""")
    parser.add_argument('--outputfilename', help="""Output file name.""")
    parser.add_argument('--header', nargs='+', type=str, help="""Space separated list of column names for output file""")

    args = parser.parse_args()

    filename = args.filename
    fieldSeparator = args.nameFieldSeparator
    outputfilename = args.outputfilename

    outdf = makedataframe(filename, fieldSeparator, outputfilename)
    if args.header is not None:
        outdf.to_csv(outputfilename, header=args.header, index=False, sep='\t', na_rep="NA")
    else:
        outdf.to_csv(outputfilename, header=None, index=False, sep='\t', na_rep="NA")
