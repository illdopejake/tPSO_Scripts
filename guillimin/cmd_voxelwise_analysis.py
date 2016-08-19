import sys
import argparse
import bootstrap_spatial_correlations as bsc

if __name__ == '__main__':
    """command line caller for bsc function voxelwise analysis. Run function
    without defaults for help"""

    parser = argparse.ArgumentParser(
        description='runs voxelwise correlation on a list of 3D volumes')

    ###INPUTS###

    parser.add_argument('parin',type=str,nargs=1,
        help='Path to input spreadsheet (or use -scans and -vals)')

    parser.add_argument('outfl',type=str,nargs=1,
        help = 'label for the results volume')

    parser.add_argument('outdir',type=str,nargs=1,
        help = 'path to desired output directory')


    ###OPTIONS###

    parser.add_argument('-scans',type=list,default=[],
        help='list of paths to subject scans')

    parser.add_argument('-vals',type=list,default=[],
        help='a list of values for the predictor variable')

    parser.add_argument('-out',type=str,default='r',choices=['r','t'],
        help='should script generate a tmap or rmap?')

    parser.add_argument('-nonpar',default=False,
        help='run spearmans instead of pearson correlation')

    parser.add_argument('-tid',default='',
        help='used to keep track of id for parallelization')

    parser.add_argument('-par',default=True,
        help='set up for parallelization and copy input files')

    parser.add_argument('-inter',default=False,
        help='save 4D volume used for analysis?')

    parser.add_argument('-indata',default=False,
        help='path to existing 4D data')

    ###RUN###

    if len(sys.argv) < 3:
        parser.print_help()
        print help(bsc.voxelwise_analysis)
    else:
        args=parser.parse_args()
        outstr=bsc.voxelwise_analysis(scans=args.scans,pv_vals=args.vals,outfl=args.outfl[0],outdir=args.outdir[0],out_tp=args.out,nonpar=args.nonpar,taskid=args.tid,parallel=args.par,parin=args.parin[0],indata=args.indata,intermed=args.inter)
        print outstr
