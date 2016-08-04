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
            scans = glob(os.path.join(outpth,'%s_*'%(cde)))

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
    outstr = outstr+'.nii'

    # clean up
    data = None
    os.system('rm %s'%(os.path.join(outpth,'%s_*'%(cde))))

    return outstr

def spatial_correlation_searchlight_from_NIAK_GLMs(indir,cov_img,outdir,contrast,templ_str,scalestr='scale',norm=False,eff='t',poly=1,taskid='',save=False):
    '''
    Given a directory containing a) NIAK glm.mat files at multiple scales, b) atlases at
    the same resolutions, and a covariate image, this function will do the
    following:
        1) Convert covariate image to resolution of each atlas
        2) Perform spatial correlations between a) average connectivity matrix
        seed from every region at every scale, b) the covariance image
        3) Output a dictionary summarizing the results, where key = resolution
        and value = a dataframe containing statistics for each test


    indir = path to directory containing NIAK glm.mat files and atlases.

    cov_img = path to covariate image

    outdir = path to desired output directory

    contrast = label of the desired contrast to extract average connectivity
    values from glm.mat file

    templ_str = search string to locate atlases

    scale_str = the string preceding the number indicating atlas resolution in
    the glm and atlas pahts

    norm = If set to a path, will normalize cov_img to target image that path
    points to. Uses flirt with nearest neighbour interpolation.  If False, no
    normalization will occur

    eff = if set to 'r', will take values within the 'eff' structure of the
    glm.mat file. If set to 't' will take values within the 'ttest' structure
    of the glm.mat file

    poly = If int > 1, results will include test statistics modeled with a
    n-order polynomial, where n = poly

    taskid = used to keep track of id in the case of multiple bootstrap samples

    save = If set to a string, will write results from each resolution to a spreadsheet
    with a file name indicated be string input
    '''

    cde = wr.codegen(6)

    if norm:
        print 'normalizing tmap to fmri space'
        os.system('flirt -interp nearestneighbour -in %s -ref %s -out %s/%s_rtmap'%(cov_img,norm,outdir,cde))
        cov_img = os.path.join(outdir,'%s_rtmap.nii.gz'%(cde))

    glmz = glob(os.path.join(indir,'glm*.mat'))
    dfz = {}
    for glm in glmz:
    ####################the 2 lines below suck######################
        scale = int(os.path.split(glm)[1].rsplit('_%s'%(scalestr))[1].rsplit('.')[0])
        scale_templ = glob(os.path.join(indir,'%s*%s.*'%(templ_str,scale)))[0]
        df, rdf, scalar = wr.get_parcelwise_correlation_map_from_glm_file(outdir,glm,scale_templ,scale,contrast,eff=eff,cov_msk='',cov_img=cov_img,conndf = '',poly=poly)
        dfz.update({scale: rdf})
        if save:
            rdf.to_excel(os.path.join(outdir,'%s_scl%s_res.xls'%(save,scale)))

    resdf = pandas.DataFrame(columns =['scale','parcel','measure','value','pvalue'])
    os.system('rm %s'%(os.path.join(outdir,'%s_*'%(cde))))

    return dfz

def id_sig_results(dfz,outdir,outfl,perc_thr=False,thr_tp='fwe',thr='',out_tp='samp',res_tp='all',master_ss=False,master_thr='',taskid=''):
    '''given a dict outputted from the searchlight funcion, and thresholding
    information, this function will save the top results from the searchlight
    into a run specific spreadsheet and/or a master spreadsheet (across
    bootstrap samples)

    dfz = a dict outputted from the searchlight function where the key is atlas
    resolution and value is a specific dataframe of results

    outdir = the desired output directory

    outfl = the name of the  output file (NOTE: taskid will be
    automatically appended to this)

    perc_thr: 'perc' = set threshold by top percentage of values (values determined by thr_tp)
              'top' = keep on the top result (value determined by res_tp)
              'fwe' = threshold based on # of comparisons (bonferonni)

    out_tp =  'samp' = only create output spreadsheet for the input
              'mast' = only write results to a master spreadsheet
              'both' = do both

    thr_tp =  'p' = set threshold by pvalue
              'r' = set threshold by test statistic
              'rabs' = set threshold by absolute value of test statistic'

    thr = desired significance threshold if thr_tp set to 'p','r', or if perc_thr set to 'perc'

    res_tp =  'r' = only report values from standard pearson r results
              'rho' = only report values from spearman r results
              'poly' = only report values from polynomial results
              'all' = report values from all three types of tests

    master_ss = path to the master spreadsheet. If no path is given, a master
    spreadsheet will be created. Will be ignored if out_tp is set to 'samp'

    master_thr = threshold for master spreadsheet results. If blank, will use
    value from thr. Will be ignored if out_to set to 'samp'. NOTE: master_thr
    MUST be equal to or more conservative than thr

    taskid = used to keep track of id in the case of multiple bootstrap samples

    WARNING: To avoid giant files and potentially breaking things, consider
    setting thr to very conservative stats, especially if res_tp = 'all'. Even
    'fwe' will often yield >5% of results significant. Using perc and low
    values, or r and a high thr, is recommended.

    NOTE: Function does not currently support searching for results smaller
    than a given statistic. i.e. if thr_tp and res_tp are both set to 'r', the
    functions can not search for results LESS THAN r, only greater. But see
    'rabs' option for thr_tp.
    '''

    # check inputs
    acc_vals = ['r','rho','poly', 'all']
    if res_tp not in acc_vals:
        raise ValueError('res_tp must be set to an appropriate value: r. rho, poly, or all. See documentation for id_sig_results for help ')
    acc_vals = ['samp','mast','both']
    if out_tp not in acc_vals:
        raise ValueError('out_tp must be set to an appropriate value: samp, mast, or both. See documentation for id_sig_results for help')
    acc_vals = ['p','r','rabs']
    if thr_tp not in acc_vals:
        raise ValueError('res_tp must be set to an appropriate value: p, r, or rabs. See documentation for id_sig_results for help')
    if perc_thr:
        acc_vals = ['perc','top','fwe']
        if perc_thr not in acc_vals:
            raise ValueError('perc_thr must be set to an appropriate value: perc,top,or fwe, or should be set to False. See documentation for id_sig_results for help')
    # wrangle spreadsheets

    resdf = pandas.DataFrame(columns =['scale','parcel','measure','value','pvalue'])

    if out_tp != 'samp':
        if os.path.isfile(master_ss):
            if master_ss[-3:] == 'xls' or master_ss[-4:] == 'xlsx':
                mdf = pandas.ExcelFile(master_ss).parse('Sheet1')
            elif master_ss[-3:] == 'csv':
                mdf = pandas.read_csv(master_ss)
            else:
                raise IOError('input spreadsheet filetype not recognized')
        else:
            master_ss = os.path.join(outdir,'master_ss.xls')
            print 'no master_ss found. Creating new master_ss at %s'%(master_ss)
            mdf = pandas.DataFrame()

    # determine thresholds    
    if out_tp != 'samp':
        if master_thr == '':
            master_thr = thr

        if thr_tp == 'r' or thr_tp == 'rabs':
            if master_thr < thr:
                raise ValueError('master_thr must be more conservative than thr')
        else:
            if master_thr > thr:
                raise ValueError('master_thr must be more conservative than thr')

    res_dict = {'r': 0, 'rho': 2, 'poly': 4}

    if perc_thr:
        comps = 0
        for scl,scf in dfz.iteritems():
            comps = comps + len(scf)
        if res_tp = 'all':
            comps = comps * 3

        if perc_thr == 'fwe':
            thr = (0.05/comps)
            print '%s total comparisons. Setting pvalue threshold to %s'%(comps,thr)
        else:
            vec = []
            for scl,scf in dfz.iteritems():
                for x,y in scf.iterrows():
                    if thr_tp == 'p'
                        if res_tp != 'all': 
                            vec.append(y[res_dict['res_tp']+1])
                        else:
                            for k,v in res_dict.iteritems():
                                vec.append(y[v+1])
                    elif thr_tp == 'r' or thr_tp == 'rabs':
                        if res_tp != 'all':
                            vec.append(y[res_dict['res_tp'])
                        else:
                            for k,v in res_dict.iteritems():
                                vec.append(y[v])

            if perc_thr == 'top':
                print 'acquiring top result'
                if res_tp == 'p':
                    thr = sorted(vec)[0]
                elif res_tp == 'r' or res_tp == 'rabs':
                    thr = sorted(vec)[-1]
            else:
                frac = int(comps * thr)
                if res_tp == 'p':
                    thr = sorted(vec)[frac]
                    print 'acquiring top results, r > %s'%(thr)
                else:
                    thr = sorted(vec)[-(frac+1)]
                    print 'acquiring top results, p < %s'%(thr)
    res_count = 0
    for scl,scf in dfz.iteritems():
        for x,y in sdf.iterrows():
            if thr_tp == 'p'
                if res_tp != 'all':
                    if y[(res_dict[res_tp]) + 1] < thr:
                        update_spreadsheet(resdf,scl,x,y,res_tp)
                        if out_tp != 'samp':
                            if master_thr == thr:
                                update_spreadsheet(resdf,scl,x,y,res_tp,v='',df_tp='mast',res_count=res_count,taskid=taskid)
                            else:
                                if y[(res_dict[res_tp]) + 1] < master_thr
                                    update_spreadsheet(resdf,scl,x,y,res_tp,v='',df_tp='mast',res_count=res_count,taskid=taskid)
                else:
                    for k,v in res_dict.iteritems():
                        if y[v+1] < thr:
                            update_spreadsheet(resdf,scl,x,y,res_tp,v=v)
                                if out_tp != 'samp':
                                    if master_thr == thr:
                                        update_spreadsheet(resdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
                                    else:
                                        if y[v+1] < master_thr:
                                            update_spreadsheet(resdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
            else:
                if res_tp != 'all':
                    if thr_tp == 'rabs':
                        if abs(y[(res_dict[res_tp])) > thr:
                            update_spreadsheet(resdf,scl,x,y,res_tp)
                            if out_tp != 'samp':
                                if master_thr == thr:
                                    update_spreadsheet(resdf,scl,x,y,res_tp,v='',df_tp='mast'.res_count=res_count,taskid=taskid)
                                else:
                                    if abs(y[(res_dict[res_tp])) > thr:
                                        update_spreadsheet(resdf,scl,x,y,res_tp,v='',df_tp='mast'.res_count=res_count,taskid=taskid)
                    else:
                        if y[(res_dict[res_tp]) > thr:
                            update_spreadsheet(resdf,scl,x,y,res_tp)
                            if out_tp != 'samp':
                                if master_thr == thr:
                                    update_spreadsheet(resdf,scl,x,y,res_tp,v='',df_tp='mast'.res_count=res_count,taskid=taskid)
                                else:
                                    if y[(res_dict[res_tp]) > thr:
                                        update_spreadsheet(resdf,scl,x,y,res_tp,v='',df_tp='mast'.res_count=res_count,taskid=taskid)
                else:
                    for k,v in res_dict.iteritems():
                        if thr_tp == 'rabs':
                            if abs(y[v+1]) > thr:
                                update_spreadsheet(resdf,scl,x,y,res_tp,v=v)
                                if out_tp != 'samp':
                                    if master_thr == thr:
                                        update_spreadsheet(resdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
                                    else:
                                        if abs(y[v+1]) > thr:
                                            update_spreadsheet(resdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
                        else:
                            if y[v+1] > thr:
                                update_spreadsheet(resdf,scl,x,y,res_tp,v=v)
                                if out_tp != 'samp':
                                    if master_thr == thr:
                                        update_spreadsheet(resdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
                                    else:
                                        if y[v+1] > thr:
                                            update_spreadsheet(resdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)

    if out_tp != 'mast':
        outstr = os.path.join(outdir,'%s%s.xls'%(outfl,taskid))
        print 'writing results to spreadsheet at %s'%(outstr)
        resdf.to_excel(outstr)

    if out_tp != 'samp':
        print 'writing results to master spreadsheet at %s'%(master_ss)
        mdf.to_excel(master_ss)

def update_spreadsheet(indf,scl,x,y,res_tp,v='',df_tp='samp',res_count=0,taskid=''):
    '''update spreadsheet with values indexed specifically from searchlight
    generated spreadsheet'''

    res_dict = {'r': 0, 'rho': 2, 'poly': 4}

    if df_tp == 'samp':
        if res_tp == 'all':
            indf.ix['%s: %s'%(scl,x), 'value'] = y[v]
            indf.ix['%s: %s'%(scl,x), 'pvalue'] = y[v+1]
        else:
            indf.ix['%s: %s'%(scl,x), 'value'] = y[res_dict[res_tp]]
            indf.ix['%s: %s'%(scl,x), 'pvalue'] = y[res_dict[res_tp]+1]
        indf.ix['%s: %s'%(scl,x), 'scale'] = scl
        indf.ix['%s: %s'%(scl,x), 'parcel'] = x
        indf.ix['%s: %s'%(scl,x), 'measure'] = res_tp


    elif df_tp == 'mast':
        indf.ix['samp_%s'%(taskid),'res%s_scale'%(res_count)] = scl
        indf.ix['samp_%s'%(taskid),'res%s_parcel'%(res_count)] = x
        indf.ix['samp_%s'%(taskid),'res%s_measure'%(res_count)] = res_tp
        if res_tp == 'all':
            indf.ix['samp_%s'%(taskid),'res%s_value'%(res_count)] = y[v]
            indf.ix['samp_%s'%(taskid),'res%s_pvalue'%(res_count)] =y[v+1]
        else:
            indf.ix['samp_%s'%(taskid),'res%s_value'%(res_count)] = y[res_dict[res_tp]]
            indf.ix['samp_%s'%(taskid),'res%s_pvalue'%(res_count)] = y[res_dict[res_tp]+1]
        res_count = res_count + 1

        return res_count
#def replicate_previous_findings
