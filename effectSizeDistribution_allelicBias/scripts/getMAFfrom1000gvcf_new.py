import subprocess as sp
import os
import argparse
import glob
import sys
import shutil
import pandas
import csv
import errno

def printWait(workflowFile):
    workflowFile.write("\n# drmr:wait\n")
    
def printInfo(workflowFile, info):
    workflowFile.write("\n###\n # %s \n###\n" % (info))
    
def printResources(workflowFile, nodes, cores, memory, time):
    workflowFile.write("\n# drmr:job nodes=%s processors=%s processor_memory=%s time_limit=%s\n" %(nodes, cores, memory, time))
    
def withIonIce(x):
    return "ionice -c2 -n7 " + x

def newMkdir(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise
        
def printBash(x):
    return x.replace("\t", "\\t")

def getMaf(workflowFile, intermediateFileDir, vcfDirPath, population, inputfile, header, outputfile):
    filedf = pandas.read_csv(inputfile, sep='\t', header=0, names=header)
    
    for name, group in filedf.groupby('chrom'):
        filename = os.path.join(intermediateFileDir, "%s.%s" %(name, os.path.basename(inputfile)))
        vcfFileName = os.path.join(vcfDirPath, "ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" %(name))
        # vcfFileName = os.path.join(vcfDirPath, "ALL.%s.NoDup.recode.vcf.gz" %(name))
        df = group[['chrom','start','end']].drop_duplicates()
        df.loc[:,'chrom'] = df['chrom'].str.replace("chr","")
        df.to_csv(filename, sep='\t', header=False, index=False)
        getInfo = "--get-INFO %s " %(" --get-INFO ".join(population))
        cmd_vcftools = """vcftools --gzvcf %s --bed %s %s --out %s""" %(vcfFileName, filename, getInfo, filename)
        cmd_chrEdit = """awk '{print "chr"$0}' OFS='\t' %s.INFO | grep -v chrCHROM > temp; mv temp %s.INFO""" %(filename, filename)
        cmd_mergeEqtlDf = """python ~arushiv/toolScripts/merge2FilesByColnames.py %s %s.INFO %s.eqtlout -h1 %s -h2 chrom end ref_x alt altFreq -on chrom end -t inner \n""" %(inputfile, filename, filename,  " ".join(header))
        workflowFile.write(printBash("%s ; %s ; %s \n" %(cmd_vcftools, cmd_chrEdit, cmd_mergeEqtlDf)))

    printWait(workflowFile)
    printResources(workflowFile, 1, 1, 4000, "15:00")
    string = os.path.join(intermediateFileDir, "*%s.eqtlout" %(os.path.basename(inputfile)))
    cmd_mergeOutputs = """cat %s | sort -r | uniq  > %s """ %(string, outputfile)
    workflowFile.write(printBash(cmd_mergeOutputs))
        

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Workflow for generating drmr file to get MAF from 1000g vcfs', usage='python getMAFfrom1000gvcf.py inputfile outputfile ')
    parser.add_argument('inputfile', type=str, help=""" Input bed file. Default considers no header and bed format. """)
    parser.add_argument('outputfile', type=str, help=""" output bed file which will have ref, refFreq, alt, altFreq columns appended.""")
    parser.add_argument('-pop','--population', type=str, nargs='+', default=['EUR_AF'], help="""Population to get allele frequency for. Options: 'AF'=calulate for all samples in vcf. 'EAS_AF' = East Asian; 'AMR_AF' = Ad Mixed American; 'EUR_AF' = European; 'SAS_AF' = South Asian' (default: "EUR_AF"). Can provide multiple options""")
    parser.add_argument('-h1', '--file1header', nargs='+', default = ['chrom','start','end'], help = """Provide header names for bed file. Requires first 3 columns named = ['chrom','start','end']""")
    parser.add_argument('-name','--workflowFile', type=str, nargs='?', default="workflow", help="""Name string for current workflow. Will be appended to files generated (default: "workflow")""")
    parser.add_argument('--now', action='store_const', const='now', default='wait', help="""Submit drmr job file now.""")
    parser.add_argument('-vcf', '--vcfDirPath', type=str, default='/lab/data/genomes/human/hg19/1000GenomesDownloads/', help="""Path to directory containing duplicates removed vcf files. Default = /lab/data/genomes/human/hg19/1000GenomesDownloads/""")
    parser.add_argument('-od', '--outputFileDir', type=str, default='intermediateFiles', help="""Directory for intermediate output files.""")
    args = parser.parse_args()
    
    workflowFile = args.workflowFile
    vcfDirPath = args.vcfDirPath
    population = args.population
    inputfile = args.inputfile
    header = args.file1header
    outputfile = args.outputfile
    intermediateFileDir = args.outputFileDir

    
    with open(workflowFile, 'w') as f:
        
        printInfo(f, "Get MAF using 1000g VCF files")
        printResources(f, 1, 1, 4000, "45:00")
        getMaf(f, intermediateFileDir, vcfDirPath, population, inputfile, header, outputfile)
        

        
    if args.now != "wait":
        sp.call("drmr %s" % workflowFile, shell=True)
