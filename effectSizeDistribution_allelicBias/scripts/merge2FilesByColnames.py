import pandas
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge 2 dataframes by column names', usage='python merge2FilesNyColnames.py file1.dat file2.dat outputfile.dat -c1 col1 col2 col3 -c2 head1 head2 head3')
    parser.add_argument('file1', help="""File 1, with header names""")
    parser.add_argument('file2', help="""File 2, with header names""")
    parser.add_argument('outputfile', help ="""Output file.""")
    parser.add_argument('-h1', '--file1header', nargs='+', help = """If file1 doesn't have headers, provide column names""")
    parser.add_argument('-h2', '--file2header', nargs='+', help = """If file2 doesn't have headers, provide column names""")
    parser.add_argument('-on', '--mergeon', nargs='+', help = """Column names to merge on if they match between files""")
    parser.add_argument('-c1', '--file1colnames', nargs='+', help = """Column headers to merge on for file1""")
    parser.add_argument('-c2', '--file2colnames', nargs='+', help = """Column headers to merge on for file2""")
    parser.add_argument('-sep1', '--separatorfile1', default="\t", help = """File 1 field separator. Default = tab""")
    parser.add_argument('-sep2', '--separatorfile2', default="\t", help = """File 2 field separator. Default = tab""")
    parser.add_argument('-t', '--mergetype', default="inner", help = """left, right, outer or inner join. Default='inner'""")
    parser.add_argument('-k', '--keepCols', nargs='+', help = """Keep these columns in the final output. If some columns not provided in -c1 or -c2 have same names in two files, remember to include _x or _y""")
    parser.add_argument('-rdup', '--removeDuplicates', action = 'store_true', help = """Remove duplicates considering all columns""")
    args = parser.parse_args()

    if args.file1header is not None:
        d1 = pandas.read_csv(args.file1, sep = args.separatorfile1, header = None, names = args.file1header)
    else:
        d1 = pandas.read_csv(args.file1, sep = args.separatorfile1)

    if args.file2header is not None:
        d2 = pandas.read_csv(args.file2, sep = args.separatorfile2, header = None, names = args.file2header)
    else:
        d2 = pandas.read_csv(args.file2, sep = args.separatorfile2)
    
    col1 = args.file1colnames
    col2 = args.file2colnames
    mergetype=args.mergetype
    print(d1.columns)
    print(d2.columns)

    if args.mergeon is not None:
        d = pandas.merge(d1, d2, how=mergetype, on=args.mergeon)
    else:
        d = pandas.merge(d1, d2, how=mergetype, left_on=col1, right_on=col2) 

    if args.keepCols is not None:
        d = d[args.keepCols]

    if args.removeDuplicates:
        d.drop_duplicates(inplace=True)
        
    d.to_csv(args.outputfile, sep='\t', index=False, na_rep="NA")
