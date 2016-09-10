import sys
import argparse
import bootstrap_spatial_correlations as bsc


if __name__ == '__main__':
    """command line caller of the bsc function collect results. Run
    function without defaults for help"""

    parser = argparse.ArgumentParser(
        description='collects summary statistics into a spreadsheet')

    ###INPUTS###
    parser.add_argument('ssdir',type=str,nargs=1,
        help='the directory where the results spreadsheets are')

    parser.add_argument('ssstr',type=str,nargs=1,
        help='unique string identifying results spreadsheets.')

    parser.add_argument('ssext',type=str,nargs=1,
        help='string representing extension of results spreadsheet (e.g. csv)')

    parser.add_argument('outdir',type=str,nargs=1,
        help='path for output files')

    ###OPTIONS###

    parser.add_argument('-outfl',type=str,default='',
        help='string label to identify output spreadsheet, if desired')

    parser.add_argument('-thr',type=str,default='r',choices=['r','rabs','p'],
        help='collect results based on what statistical value')

    parser.add_argument('-resamp',default=False,
        help='If True, 0<float<1 representing perc. of value distrib. 2 report')

    parser.add_argument('-test',default='',
        help='Dataframe or path to spreadsheet used to test bootstrap results')

    parser.add_argument('-sum',default=False,
        help = 'If True, will output somewhat pointless summary statistics')

    ###RUN###

    if len(sys.argv) < 4:
        parser.print_help()
        print help(bsc.collect_results)
    else:
        args = parser.parse_args()
        bsc.collect_results(ss_dir=args.ssdir[0],ss_str=args.ssstr[0],ss_ext=args.ssext[0],outdir=args.outdir[0],outfl=args.outfl,thr_tp=args.thr,resample=args.resamp,permtest=args.test,summary = args.sum)

