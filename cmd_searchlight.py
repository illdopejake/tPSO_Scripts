import sys
import argparse
import bootstrap_spatial_correlations as bsc

if __name__ == '__main__':
    """run spatial correlations between a 3D volume and every connectivity map
    at mutiple resolutions"""

    parser = argparse.ArgumentParser(
            description='runs multi-input spatial correlations')

    ###INPUTS###

    parser.add_argument('indir',type=str,nargs=1,
        help='Path to input glm.mat files and and atlases')

    parser.add_argument('cov_img',type=str,nargs=1,
        help='Path to covariate image')

    parser.add_argument('outfl',type=str,nargs=1,
        help='str to label results spreadsheets')

    parser.add_argument('outdir',type=str,nargs=1,
        help='Desired output directory')

    parser.add_argument('contrast',type=str,nargs=1,
        help='label for glm file structure housing connectivity values')

    parser.add_argument('templ_str',type=str,nargs=1,
        help='search string to locate atlases')

    ###OPTIONS###

    parser.add_argument('-sclstr',type=str,default='scale',
        help='atlas family name (str proceeding resolution number')

    parser.add_argument('-norm',default=False,
        help='If path, norm cov_img to fMRI space using path as template')

    parser.add_argument('-eff',type=str,default='r',choices=['r','t'],
        help='extract values from eff or ttest structure in glm file?')

    parser.add_argument('-poly',type=int,default=1,
        help= 'where n=poly, model with n-ordered polynomial')

    parser.add_argument('-tid',default='',
        help='used to keep track of id for parallelization')

    ###RUN###

    if len(sys.argv) < 6:
        parser.print_help()
        print help(bsc.spatial_correlation_searchlight_from_NIAK_GLMs)
    else:
        args=parser.parse_args()
        dfz = bsc.spatial_correlation_searchlight_from_NIAK_GLMs(indir=args.indir[0],cov_img=args.cov_img[0],outdir=args.outdir[0],contrast=args.contrast[0],templ_str=args.templ_str[0],scalestr=args.sclstr,norm=args.norm,eff=args.eff,poly=args.poly,taskid=args.tid,save=args.outfl[0])



