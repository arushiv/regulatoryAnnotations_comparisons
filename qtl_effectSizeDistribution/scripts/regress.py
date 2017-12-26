import pandas
import numpy
import statsmodels.formula.api as sm
import argparse

# chrom   proxyStart      proxyEnd        indexPos        gene_id gene_name       tss_distance    ref     alt     snp     maf     pval_nominal    slope   qval    Position        Alleles cell    annotation
def getOpts():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputfile', help="""name of the inputfile""")
    parser.add_argument('--outputfile', help="""name of the outputfile""")
    parser.add_argument('--typereg', help="""choose from [allpeaks, intpeak]""")

    return parser
    

""" best association per SNP """
def bestAssoc(d):
    d1 = d.sort_values('BetaAdjustedPval')
    d1 = d.drop_duplicates(['cs','es','is'])
    return d1

"""remove duplicates"""

def statevar_regress(d):
    d1 = d.drop_duplicates(['chrom','indexPos','gene_name','annotation'])
    res1 = sm.ols(formula="absoluteBeta ~ absoluteDist + statevar", data=d1).fit()
    print(res1.summary())

""" dist + state + num peak overlapping LD """
def statevar_ld_regress(d):
    d1 = d.groupby(['chrom','indexPos','gene_name','annotation', 'absoluteBeta', 'absoluteDist', 'statevar']).size().reset_index(name='numLD')
    # print(d1)
    res1 = sm.ols(formula="absoluteBeta ~ absoluteDist + statevar + numLD", data=d1).fit()
    print(res1.summary())

""" dist + state + num total LD """
def statevar_totalld_regress(d):
    res1 = sm.ols(formula="absoluteBeta ~ absoluteDist + statevar + totalnumLD", data=d).fit()
    print(res1.summary())

""" regress with atac peak length """

def peakLength_regress(d):

    # d1 = d.drop_duplicates(['cs','is','GeneName','chromatinState'])
    d1 = d.drop_duplicates(['cs','is','GeneName','chromatinState','peakChr','peakStart','peakEnd'])
    res1 = sm.ols(formula="absoluteBeta ~ absoluteDist + peakLength", data=d1).fit()
    print(res1.summary())

def fixForChromStates(d):
    d = d[d['annotation'].isin(['stretchEnhancer','hotRegions'])]
    # d = d[d['annotation'].isin(['stretchEnhancer','typicalEnhancer'])]

    d.loc[:,'statevar'] = numpy.where(d['annotation'] == "stretchEnhancer", 1, 0)
    return d

def forAllContigEnhancers(d):
    d1 = d.drop_duplicates(['cs','is','GeneName','enhancerChr', 'enhancerStart', 'enhancerEnd'])
    res1 = sm.ols(formula="absoluteBeta ~ absoluteDist + intersectPeakLength ", data=d1).fit()
    print(res1.summary())

if __name__ == '__main__':
    parser = getOpts()
    args = parser.parse_args()

    d = pandas.read_csv(args.inputfile, sep='\t')

    d.loc[:,'absoluteBeta'] = abs(d['slope'])
    d.loc[:,'absoluteDist'] = abs(d['tss_distance'])
    
    if args.typereg == "intpeaks":
        d = fixForChromStates(d)
        statevar_regress(d)
        statevar_ld_regress(d)
        statevar_totalld_regress(d)
        
    elif args.typereg == "allpeaks":
        forAllContigEnhancers(d)

    else:
        print("select valid type")
