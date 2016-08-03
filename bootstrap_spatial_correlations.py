import os
from glob import glob
import pandas
import random, math
import numpy as np
import nibabel as ni
import scipy.stats.stats as st
import propagation_correlations as wr

def define_bootstrap_sample(ss,subcol,subpth,pv,sample_perc=0.5,num_gps=3,taskid = ''):
    """takes an existing sample of subjects from a spreadsheet and creates a
    subsample, balanced by a variable. Outputs a subsamble membership txt file and
    a list of paths corresponding to subsample member scans)

    ss = a spreadsheet containing at least a list of subjects and values for a
    predictor variable

    subcol = string corresponding to the label of the column containing values for
    the subject IDs. Note these IDs should be consistent with existing file names
    in subpth. If index, set to 'index'

    subpth = a path pointing to the directory containing subject scans

    pv = string corresponding to the label of the column containing values for the
    predictor variable

    sample_perc = the percent of the total sample to be used in each subsample.
    Default is 0.5 for 50%. Must be float >0 and <=1

    num_gps = the number of groups used to balance subsample on predictor variable.
    Defualt is 3. Value must be int

    task_id = used to keep track of id in the case of multiple bootstrap samples

    ***IMPORTANT NOTE***
    Script has the following assumptions:
    1) scans in subpth directory have ids from subcol in them, but the filename of
    the scan does not /start/ with the id.
    2) There are no more than 1 scan per ID within subpth directory
    """

    # prep spreadsheet

    if type(num_gps) != int:
        raise IOError('numgps must be an integer')

    if sample_perc > 1 or sample_perc <=0:
        raise IOError('sample_perc must be a float between 0 and 1')

    if ss[-3:] == 'xls' or ss[-4:] == 'xlsx':
        subdf = pandas.ExcelFile(ss).parse('Sheet1')
    elif ss[-3:] == 'csv':
        subdf = pandas.read_csv(ss)
    else:
        raise IOError('input spreadsheet filetype not recognized. Please use .xls, .xlsx or .csv')

    if subcol != 'index':
        if subcol in subdf.columns.tolist():
            df.index = df[:][subcol].tolist()
        else:
            raise IOError('there is no column in the spreadsheet called %s'%(subcol))

    for sub in subdf.index.tolist():
        if len(glob(os.path.join(subpth,'*%s*'%(sub)))) < 1:
            print 'no scan was found for subject %s. Removing from analysis.'%(sub)
            subdf.drop(sub,axis = 0, inplace=True)

    if pv not in subdf.columns.tolist():
        raise IOError('there is no column in the spreadsheet called $s'%(pv))
    else:
        subdf = subdf.sort(pv)
    allsubs = subdf.index.tolist()

    # extract subsample

    subsamp_n = len(allsubs) / num_gps
    subsamp_dict = {}
    for i in range(1,(num_gps+1)):
        if i == 1:
            subsamp_dict.update({i: allsubs[:subsamp_n]})
        else:
            subsamp_dict.update({i: allsubs[(subsamp_n * (i-1)):(subsamp_n *
            (i))]})

    subsamp = []
    n_per_gp = int(subsamp_n * sample_perc)
    for gp,lst in subsamp_dict.iteritems():
        random.shuffle(lst)
        for i in range(n_per_gp):
          subsamp.append(lst[i])

   # collect outputs
    scans = []
    for sub in subsamp:
        scn = glob(os.path.join(subpth,'*%s*'%(sub)))[0]
        scans.append(scn)

    pv_vals = []
    for sub in subsamp:
        pv_vals.append(subdf.ix[sub,pv])

    # make membership record
    ndf = pandas.DataFrame(index=subsamp)
    if len(taskid) > 0:
        ndf.to_csv('%s_subsample_membership'%(taskid))
    else:
        cde = wr.codegen(6)
        ndf.to_csv('%s_subsample_membership'%(cde))

    return scans, pv_vals

#def generate_RSNs_for_subsample

def convert_r2t(r,samp_df):
    t = math.sqrt((samp_df * r**2) / (1 - r**2))
    if r < 0:
        t = (t*-1)

    return t

def voxelwise_analysis(scans,pv_vals,outfl,outdir,out_tp='r',nonpar=False,taskid='',parallel=False,indata=False,intermed=False):
    """Given a list of 3D volumes and a corresponding list of values for a
    predictor variable, run voxelwise correlations.

    scans = a list of paths corresponding to 3D nifti volumes

    pv_vals = a list of values corresponding to subject values, in the same
    order as scans

    outfl = string to be used to identifiy the outfile map (NOTE, taskid will
    be automatically appended to this string)

    outdir = the path to the directory where the outfiles will be written to

    out_tp = the type of map to be generated. 'r' will generate a voxelwise
    rmap, whereas 't' will generate a voxelwise tmap

    nonpar = Set to True to set voxelwise analysis from pearson to spearman

    indata = In case 4D data is already available in an existing variable

    taskid = used to keep track of id in the case of multiple bootstap samples

    parallel = if true, script will copy scans into working directory to allow
    for multiple concurrent processes.

    intermed = if true, script will not delete the 4D volume used to run the
    voxelwise analysis

    WARNING: As of now, this script is extremely computationally intensive for
    a local machine. Running a subsample of 135 subjects required 5 GB memory,
    used 100% of my CPU at times, and took well over 10 minutes...

    NOTE: As of now, script does not regress out confounding variables or mask
    analysis
    """

    cde = wr.codegen(6)

    if indata:
        data = indata
    else:
        if parallel:
            for ind,scan in scans:
                os.system('cp %s %s/%s_subject%s'%(scan,outpth,cde,ind))
            scans = glob(os.path.join(outpth,'%s_*'%(cde))

        if intermed:
            intfl = 'intfile%s'%(taskid)
        else:
            intfl = '%s_intfile'%(cde)

        # create 4D volume
        cmd = 'fslmerge -t %s'%(os.path.join(outdir,intfl))
        for scn in scans:
            cmd = cmd+' %s'%(scn)

        print 'creating input 4D volume...'
        os.system(cmd)

        print 'loading data'
        data = ni.load(os.path.join(outdir,intfl)).get_data()

    # run voxelwise analysis
    print 'beginning analysis...'
    x,y,z,t_dim = data.shape
    results = np.zeros((x,y,z))
    aff = ni.load(scans[0]).affine

    for xind in range(x):
        for yind in range(y):
            for zind in range(z):
                if all(data[xind][yind][zind][:]) == 0:
                    continue
                else:
                    dv_vals = data[xind][yind][zind][:]
                    if nonpar:
                        r,p = st.spearmanr(dv_vals,pv_vals)
                    else:
                        r,p = st.pearsonr(dv_vals,pv_vals)
                    if out_tp == 't':
                        r = convert_r2t(r,t_dim)
                    results[xind,yind,zind] = r
        print 'finished %s/%s job clusters'%(xind,x)

    # write image
    outstr = '%s%s'%(os.path.join(outdir,outfl),taskid)
    print 'writing image to %s'%(outstr)
    nimg = ni.Nifti1Image(results,aff)
    ni.save(nimg,outstr)


    return outstr

def spatial_correlation_searchlight

#def id_sig_results

#def replicate_previous_findings
