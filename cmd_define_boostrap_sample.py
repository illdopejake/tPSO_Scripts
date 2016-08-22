import sys
import argparse
import bootstrap_spatial_correlations as bsc


if __name__ == '__main__':
    """command line caller of the bsc function define bootstrap sample. Run
    function without defaults for help"""

    parser = argparse.ArgumentParser(
        description='extracts random, balanced subsample from spreadsheet')

    ###INPUTS###
    parser.add_argument('ss',type=str,nargs=1,
        help='a spreadsheet wih list of subjects and pv values')

    parser.add_argument('subcol',type=str,nargs=1,
        help='label of column containing values for the subject IDs.')

    parser.add_argument('subpth',type=str,nargs=1,
        help='a path pointing to the directory containing subject scans')

    parser.add_argument('pv',type=str,nargs=1,
        help='label of column containing values for the predictor variable')

    parser.add_argument('outpth',type=str,nargs=1,
        help='desired output directory')

    ###OPTIONS###

    parser.add_argument('-resamp',default=False,
        choices=[False,'False','subsamp','jackknife','permute','bootstrap'],
        help='pick which resampling method you wish to use, if any')

    parser.add_argument('-par',default=False,
        help='if True, makes compatible for command line based parallelization')

    parser.add_argument('-rand', default=False,
        help='Seed for PRNG, or leave False to randomly generate seed')

    parser.add_argument('-tid',default='',
        help='used to keep track of id for parallelization')

    parser.add_argument('-perc',type=float,default=0.5,
        help='the percent of the total sample to be used in each subsample.')

    parser.add_argument('-gps',type=int,default=3,
        help='# of groups used to balance subsample on predictor variable.')

    ###RUN###

    if len(sys.argv) < 5:
        parser.print_help()
        print help(bsc.define_bootstrap_sample)
    else:
        args = parser.parse_args()
        if args.par:
            scans,pv_vals,flpth=bsc.define_bootstrap_sample(args.ss[0],args.subcol[0],args.subpth[0],args.pv[0],args.outpth[0],args.resamp,args.par,args.rand,args.tid,args.perc,args.gps)
            print flpth
        else:
            scans,pv_vals=bsc.define_bootstrap_sample(args.ss[0],args.subcol[0],args.subpth[0],args.pv[0],args.outpth[0],args.resamp,args.par,args.rand,args.tid,args.perc,args.gps)

