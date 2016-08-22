import sys
import argparse
import bootstrap_spatial_correlations as bsc


if __name__ == '__main__':
    """command line caller of the bsc function create_output_images. Run
    function without defaults for help"""

    parser = argparse.ArgumentParser(
        description='create images based on results')

    ###INPUTS###
    parser.add_argument('resdf',nargs=1,
        help='a pandas results dataframe. Ignored if par')

    parser.add_argument('scldir',type=str,nargs=1,
        help='path to atlas directory.')

    parser.add_argument('sclstr',type=str,nargs=1,
        help='string label for atlases in scldir')

    parser.add_argument('outdir',type=str,nargs=1,
        help='desired output directory')

    parser.add_argument('outstr',type=str,nargs=1,
        help='label for output images')

    parser.add_argument('input',type=str,nargs=1,
        help='label for measure to be used to generate results (e.g r or rho')

    ###OPTIONS###

    parser.add_argument('-outtp',type=str,default='resamp',
        choices=['resamp','cross','ci','all'],
        help='pick which types of images to output')

    parser.add_argument('-par',default=False,
        help='if True, path pointing to input spreadsheet instead of resdf')

    ###RUN###

    if len(sys.argv) < 5:
        parser.print_help()
        print help(bsc.create_output_images)
    else:
        args = parser.parse_args()
        bsc.create_output_images(resdf=args.resdf[0],scldir=args.scldir[0],sclstr=args.sclstr[0],outdir=args.outdir[0],outstr=args.outstr[0],input_tp=args.input[0],output_tp=args.outtp,par=args.par)
            