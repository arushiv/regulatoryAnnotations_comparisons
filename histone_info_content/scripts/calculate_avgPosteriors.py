import os
import pandas
import numpy
import pybedtools
import math
import argparse
import subprocess as sp

def log_with_nan(x, y):
    try:
        return math.log(x, y)
    except ValueError:
        "Places with log(0)"
        return float('nan')
        

def calculate_avgPosteriors(out, **kwargs):
    """Concat sum od probabilities on intersected index. Then calculate average."""
    
    def fixdf(filename):
        name = os.path.basename(filename).replace(kwargs['replaceString'], "")
        d = pandas.read_csv(filename, sep='\t', header=1, usecols = kwargs['states'])
        prob_sums = d.loc[indices, :].apply(sum, axis=1)
        prob_sums.rename(name, inplace=True)
        return prob_sums
    
    outdf = pandas.concat([out] + [fixdf(filename) for filename in kwargs['inputfiles']] , axis=1, join="inner")
    avg_posteriors = outdf.groupby(['chrom','start','end']).apply(numpy.mean).drop(['start','end'], axis=1)

    # avg_posteriors.rename(columns={kwargs['cellName'] : 'avg_posterior'}, inplace=True)
    avg_posteriors.reset_index(inplace = True)

    return avg_posteriors
    

def entropy_function(x):
    y = -x * log_with_nan(x, 2)
    return y

def calculate_information(d):
    """
    Relative RPKM = x(g,t)/sum(x(g,t))
    I = Relative RPKM * (Entropy max - Entropy(g) )
    Entropy max (all probabilities equal)= sum(-1/N * log2(-/N)) = log2(N) 
    """

    N = len(d.columns)
    print(f"Number of features = {N}")
    
    entropy_max = math.log(N, 2)
    print(f"max entropy = {entropy_max}")
    
    d_relative_rpkm = d.div(d.sum(axis=1), axis=0)
    # print("relative RPKM")
    # print(d_relative_rpkm)
    d_relative_rpkm.loc[:,'entropy'] = d_relative_rpkm.applymap(lambda x: entropy_function(x)).sum(axis=1)
    # print("with entropy")
    # print(d_relative_rpkm)
    d_information = d_relative_rpkm.apply(lambda x: x*(entropy_max - x['entropy']), axis=1)
    d_information.drop('entropy', axis=1, inplace=True)
    d_information.loc[:,'change_in_entropy'] = entropy_max - d_relative_rpkm['entropy']
    # print("information")
    # print(d_information)
    return d_information


def work_on_calculating_infoContent(avg_posteriors, **kwargs):
    avg_posteriors.loc[:,'avg_post'] = avg_posteriors[kwargs['cellName']]
    avg_posteriors.set_index(['chrom', 'start', 'end', 'avg_post'], inplace=True)
        
    """
    Calculate information
    """
    out = calculate_information(avg_posteriors)
    out.reset_index(inplace=True)
    out.loc[:,'information'] = out[kwargs['cellName']]
    out.loc[:,'cell'] = kwargs['cellName']
    out.loc[:,'annotation'] = kwargs['annotName']
    
    return out


def getOpts():
    parser = argparse.ArgumentParser(description='Calculate average posterior probabilities for a feature - averaged across sum of probabilities of given chromhmm states per 200 bp tiles.'
                                     'Also option to calculate information content', usage='python calculate_avgPosteriors.py ')
    parser.add_argument('--inputfiles', required=True, nargs='+',
                        help="""Input probabilty files '_posterior.txt' from ChromHMM. Multiple files specified as space separated lists.""")
    parser.add_argument('--annotationfile', required=True,
                        help="""Annotation bed file - avg probabilities will be returned for each feature in this bed file.""")
    parser.add_argument('--tilefile', required=True,
                        help="""200 bp tile file""")
    parser.add_argument('--outputfile', required=True,
                        help ="""Output file.""")
    parser.add_argument('--states', nargs='+', required=True,
                        help = """States to sum up probabilitied over each 200 bp tile. Give a space separated list.""")
    parser.add_argument('--replaceString', required=True,
                        help = """String to replace from inputfile basename to get cell type name""")
    parser.add_argument('--annotName', required=True,
                        help = """Annotation name""")
    parser.add_argument('--cellName', required=True,
                        help = """Annotation cell type name""")
    parser.add_argument('--infoContent_file',
                        help = """Outputfile which will contain information content for cell type and annotation of the annotation file""")
    
    return parser

if __name__ == '__main__':
    parser = getOpts()
    args = parser.parse_args()

                                                                                            
    # get indices
    bedchrom = pybedtools.BedTool(args.tilefile)
    dsegment = pybedtools.BedTool(args.annotationfile)
    
    out = bedchrom.intersect(dsegment, wa=True, wb=True).to_dataframe(names=['c','s','e','index','chrom','start','end'])[['index', 'chrom', 'start', 'end']]

    indices = out['index'].tolist()
    out.set_index('index', inplace=True)
   
    if not indices:
        sp.call(f"""                                                                                                                                                                                     
        touch {args.outputfile}          
        """, shell=True)
        if args.infoContent_file is not None:
            sp.call(f"""                                                                                                                                                                       
            touch {args.infoContent_file}          
            """, shell=True)
            
    else:
        avg_posteriors = calculate_avgPosteriors(out, **vars(args))
        avg_posteriors.to_csv(args.outputfile, sep='\t', index=False)
        if args.infoContent_file is not None:
            d_infoContent = work_on_calculating_infoContent(avg_posteriors, **vars(args))
            d_infoContent.to_csv(args.infoContent_file, sep='\t', index=False, na_rep="NA")
        # avg_posteriors[['chrom','start','end','cell','annotation','avg_posterior','infoContent']].to_csv(output.segmented_posteriors, sep='\t', index=False)                             

