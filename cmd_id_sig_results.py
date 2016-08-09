ulmport sys
import argparse
import bootstrap_spatial_correlations as bsc

if __name__ == '__main__':
    """extract spatial correlation results using specified thresholds"""

    parser = argparse.ArgumentParser(
        description='extract results from spatial correlation analysis')

    ###INPUTS###

    parser.add_argument('in_str',type=str,nargs=1,
        help='string used to ID results spreadsheets')

    parser.add_argument('outdir',type=str,nargs=1,
        help='desired output directory')

   ###OPTIONS###

    parser.add_argument('-dfz',type=dict,default={},
        help='Dict where key=resolution and value=results DataFrame')

    parser.add_argument('-outfl',type=str,default='outfl',
        help='name for output file')

    parser.add_argument('-perc',default='top',choices=['top','fwe','perc',False],
        help='If true, threshold based on percentage')

    parser.add_argument('-tp',type=str,default='r',choices=['p','r','rabs'],
        help='Threshold results based on p-,r-,or absolute r-values?')

    parser.add_argument('-thr',default='',
        help='Desired statistical or percentage threshold')

    parser.add_argument('-outtp',type=str,default='samp',choices=['samp','mast','both'],
        help='Write output to sample-specific or master spreadsheet, or both?')

    parser.add_argument('-res',type=str,default='all',choices=['r','rho','poly','all'],
        help='If available, which type of value to threshold on')

    parser.add_argument('-mss',type=str,default='',
        help='Path to master spreadsheet if -outtp not set to samp')

    parser.add_argument('-mthr',default='',
        help='Thr for master spreadsheet, if different from thr')

    parser.add_argument('-save',type=bool,default=False,
        help='If True, will not delete input results spreadsheets')

    parser.add_argument('-tid',default='',
        help='used to keep track of id for parallelization')

    ###RUN###

    if len(sys.argv) < 2:
        parser.print_help()
        print help(bsc.id_sig_results)
    else:
        args=parser.parse_args()
        bsc.id_sig_results(dfz=args.dfz,outdir=args.outdir[0],outfl=args.outfl,perc_thr=args.perc,thr_tp=args.tp,thr=args.thr,out_tp=args.outtp,res_tp=args.res,master_ss=args.mss,master_thr=args.mthr,par=args.in_str[0],parsave=args.save,taskid=args.tid)


