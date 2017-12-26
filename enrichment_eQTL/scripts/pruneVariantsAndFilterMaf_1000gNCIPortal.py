
import argparse
import errno
import os
import subprocess as sp
import sys
import time
import glob

from bs4 import BeautifulSoup, SoupStrainer
import pandas

"""
Script to prune SNPs based on LD in specific populations using 1000g phase 3 data using the NCI portal https://analysistools.nci.nih.gov/LDlink/?tab=snpclip. Can also provide a MAF threshold to filter variants.
"""
    

def newmkdir(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise


def pruneNCI(d, resultFileDir, population, r2Threshold, mafThreshold, sortName ):
    if sortName is not None:
        d.sort_values(by=sortName, ascending=True, inplace=True)

    parameters = {
        'population' : population,
        'r2Threshold' : r2Threshold,
        'mafThreshold' : mafThreshold,
    }

    def printCommands(group):
        variantList = "\\n".join(group['snp'].tolist())
        new = {
            'variantList' : variantList,
            'result' : os.path.join(resultFileDir, "%s_pruningResult.dat" %(group.name))
        }

        parameters.update(new)
            
        cmd = """curl -k -H "Content-Type: application/json" -X POST -d '{{"snps": "{variantList}", "pop": "{population}", "r2_threshold": "{r2Threshold}", "maf_threshold": "{mafThreshold}"}}' 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/snpclip' > {result}""".format(**parameters)
        sp.call(cmd, shell=True)

    grouped = d.groupby('chrom')
    grouped.apply(printCommands)

def mergeResults(resultFileDir, d, outputfile):
    def getPrunedSnps(resultFile):
        print(resultFile)

        d = pandas.read_csv(resultFile, sep='\t', names=['snp', 'Position', 'Alleles', 'Details', 'extra'])
        d = d[d['Details'] == "Variant kept."]
        d = d[['snp', 'Position', 'Alleles']]
    
        return(d)


    resultdf = pandas.concat([getPrunedSnps(resultFile) for resultFile in glob.glob(os.path.join(resultFileDir, "*.dat"))])
    dout = pandas.merge(d, resultdf, on=['snp'], how="inner")
    dout.to_csv(outputfile, index=False, sep='\t', na_rep="NA")
    # return(dout)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to prune SNPs based on LD in specific populations using 1000g phase 3 data using the NCI portal https://analysistools.nci.nih.gov/LDlink/?tab=snpclip. Can also provide a MAF threshold to filter variants. Should supply a sorted list of variants by p-value or other applicable criteria to retain SNPs', usage='python pruneVariantsAndFilterMaf_1000gNCIPortal.py <variantfile> <outputfile> -p CEU -r2 0.8 -maf 0.2 ')
    parser.add_argument('variantfile', help="""file with rsIDs to search proxies of. Also needs chrom information. Tab separated. Header name of the rsID and chrom column should be 'snp', 'chrom'. """)
    parser.add_argument('-hv', '--variantfileHeader', nargs='+', help="""If variant file does not contain header, supply space separated list here. Header name of the rsID and chrom column should be 'snp', 'chrom'""")
    parser.add_argument('-p', '--population', default='CEU', help="""Population to search proxies in. Options are: \n
(AFR) African:[[
    (YRI) Yoruba in Ibadan, Nigera
    (LWK) Luhya in Webuye, Kenya
    (GWD) Gambian in Western Gambia
    (MSL) Mende in Sierra Leone
    (ESN) Esan in Nigera
    (ASW) Americans of African Ancestry in SW USA
    (ACB) African Carribbeans in Barbados ]]
(AMR) Ad Mixed American [[
    (MXL) Mexican Ancestry from Los Angeles, USA
    (PUR) Puerto Ricans from Puerto Rico
    (CLM) Colombians from Medellin, Colombia
    (PEL) Peruvians from Lima, Peru ]]
(EAS) East Asian [[
    (CHB) Han Chinese in Bejing, China
    (JPT) Japanese in Tokyo, Japan
    (CHS) Southern Han Chinese
    (CDX) Chinese Dai in Xishuangbanna, China
    (KHV) Kinh in Ho Chi Minh City, Vietnam ]]
(EUR) European [[
    (CEU) Utah Residents from North and West Europe
    (TSI) Toscani in Italia
    (FIN) Finnish in Finland
    (GBR) British in England and Scotland
    (IBS) Iberian population in Spain ]]
(SAS) South Asian [[
    (GIH) Gujarati Indian from Houston, Texas
    (PJL) Punjabi from Lahore, Pakistan
    (BEB) Bengali from Bangladesh
    (STU) Sri Lankan Tamil from the UK
    (ITU) Indian Telugu from the UK \n ]]
    provide sub or superpopulation code(s) WITHOUT parantheses.""" )
    parser.add_argument('-dir', '--resultFileDir', default='results_ldPrune', help="""Directory that will contain all raw LD prune results.""")
    parser.add_argument('-r2', '--r2Threshold', type=float, default=0.8, help="""Retain variants with this r2 or lower and paste into the output file. Default = 0.8""")
    parser.add_argument('-maf', '--mafThreshold', type=float, default=0.0, help="""Retain variants with this r2 or lower and paste into the output file. Default = 0.0""")
    parser.add_argument('-s', '--sortName', type=str, help="""Give column name to sort rsIDs by before pruning. Variants further down in the list that are in high LD with top variants will be removed""")

    parser.add_argument('outputfile', help="""name of the outputfile""")
    
    args = parser.parse_args()
    resultFileDir = args.resultFileDir
    newmkdir(args.resultFileDir)
    variantfile = args.variantfile
    population = args.population
    r2Threshold = args.r2Threshold
    mafThreshold = args.mafThreshold
    sortName = args.sortName
    outputfile = args.outputfile
    
    if args.variantfileHeader is not None:
        d = pandas.read_csv(variantfile, sep='\t', header=None, names=args.variantfileHeader)
    else:
        d = pandas.read_csv(variantfile, sep='\t')


    pruneNCI(d, resultFileDir, population, r2Threshold, mafThreshold, sortName)
    mergeResults(resultFileDir, d, outputfile)
