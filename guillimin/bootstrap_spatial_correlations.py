#!/sb/home/jvogel/local_python/
import os
from glob import glob
import pandas
import random, math
import numpy as np
import nibabel as ni
import scipy.stats.stats as st
import propagation_correlations as wr
import randomstate.prng.mrg32k3a as rnd

def define_bootstrap_sample(ss,subcol,subpth,pv,outpth,sample_perc=0.5,num_gps=3,par=False,rand=False,taskid = ''):
    """takes an existing sample of subjects from a spreadsheet and creates a
    subsample, balanced by a variable. Outputs a subsamble membership txt file and
    a list of paths corresponding to subsample member scans)

    ss = a spreadsheet containing at least a list of subjects and values for a
    predictor variable

    subcol = string corresponding to the label of the column containing values for
    the subject IDs. Note these IDs should be consistent with existing file names
    in subpth. If index, set to 'index'

    subpth = a path pointing to the directory containing subject scans

    outpth = desired output directory

    pv = string corresponding to the label of the column containing values for the
    predictor variable

    sample_perc = the percent of the total sample to be used in each subsample.
    Default is 0.5 for 50%. Must be float >0 and <=1

    num_gps = the number of groups used to balance subsample on predictor variable.
    Defualt is 3. Value must be int

    par = if True, makes compatible for command line based parallelization

    task_id = used to keep track of id in the case of multiple bootstrap samples

    rand = Determine how psuedorandom generator is seeded for the randomization
    of the sample. Leave as False to use random.shuffle with a random seed, 
    which will create unreproducible samples and is not recommended for 
    parallelization. Set as an int to use int as seed for mrg322ka PRNG 
    (recommended for parallelization or reproducible results)

    Outputs a list of paths and a vector containing the values for the
    predictor variable

    ***IMPORTANT NOTE***
    Script has the following assumptions:
    1) scans in subpth directory have ids from subcol in them, but the filename of
    the scan does not /start/ with the id.
    2) There are no more than 1 scan per ID within subpth directory
    """

    # prep spreadsheet

    par = check_bool(par)
    rand = check_bool(rand)

    # PUTTING THIS IN FOR NOW TILL GUILLIMIN STOPS BEING A SHIT
    if rand:
	rand = False

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
            subdf.index = subdf[:][subcol].tolist()
        else:
            raise IOError('there is no column in the spreadsheet called %s'%(subcol))

    for sub in subdf.index.tolist():
        if len(glob(os.path.join(subpth,'*%s*'%(sub)))) < 1:
            print 'no scan was found for subject %s. Removing from analysis.'%(sub)
            #this was not compatible with guillimin's version of pandas
            #subdf.drop(sub,axis = 0, inplace=True)
            subdf=subdf.drop(sub,axis=0)

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
    for j,dictees in enumerate(subsamp_dict.iteritems()):
        gp = dictees[0]
        lst = dictees[1]
        if not rand:
            random.shuffle(lst)
        else:
            if not taskid:
                lst = randomize_4_montecarlo(rand,random.randint(1,1000000),((j+1)*10000),lst)
            else:
                lst = randomize_4_montecarlo(rand,taskid,((j+1)*10000),lst)
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
        ndf.to_csv(os.path.join(outpth,'%s_subsample_membership'%(taskid)))
    else:
        cde = wr.codegen(6)
        ndf.to_csv(os.path.join(outpth,'%s_subsample_membership'%(cde)))

    if par:
        flpth = parallel_out('dbc',outpth,scans,pv_vals,taskid)
        return scans,pv_vals,flpth
    else:
        return scans, pv_vals

#def generate_RSNs_for_subsample

def randomize_4_montecarlo(seed,jumpseed,jump_mult,rsamp):

    intp = type(rsamp[0])
    rs = rnd.RandomState(seed)
    jumpseed = jumpseed * jump_mult
    rs.jump(jumpseed)
    rand_array = rs.randn(len(rsamp))

    randf = pandas.DataFrame(np.zeros((len(rsamp),2)),columns=['val','rand'])
    for i in range(len(rsamp)):
        randf.loc[randf.index[i],'val'] = rsamp[i]
        randf.loc[randf.index[i],'rand'] = rand_array[i]
    randf = randf.sort('rand')
    nrsamp = randf[:]['val'].tolist()

    if type(nrsamp[0]) != intp:
        nrsamp = [intp(x) for x in nrsamp]

def convert_r2t(r,samp_df):
    t = math.sqrt((samp_df * r**2) / (1 - r**2))
    if r < 0:
        t = (t*-1)

    return t

def voxelwise_analysis(scans,pv_vals,outfl,outdir,out_tp='r',nonpar=False,taskid='',parallel=False,parin='',indata=False,intermed=False):
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

    indata = In case 4D data is already available in an existing variable, data
    can be specified here. Or, if an int file exists, simply add the path here.

    taskid = used to keep track of id in the case of multiple bootstap samples

    parallel = if true, script will copy scans into working directory to allow
    for multiple concurrent processes. Will also make script compatible for
    command line based parallelization.

    parin = input path for parallelization

    intermed = if true, script will not delete the 4D volume used to run the
    voxelwise analysis


    Outputs a path pointing to the newly created tmap


    WARNING: As of now, this script is extremely computationally intensive for
    a local machine. Running a subsample of 135 subjects required 5 GB memory,
    used 100% of my CPU at times, and took well over 10 minutes...

    NOTE: As of now, script does not regress out confounding variables or mask
    analysis
    """

    nonpar = check_bool(nonpar)
    parallel = check_bool(parallel)
    indata = check_bool(indata)
    intermed = check_bool(intermed)

    cde = wr.codegen(6)

    if parallel:
        scans,pv_vals = parallel_in('va',outdir,parin,cde)

    if indata:
        if type(indata) != str:
            data = indata
        else:
            #if not parallel:
            print 'loading data...'
            data=ni.load(indata).get_data()
    else:
        if intermed:
            intfl = 'intfile%s'%(taskid)
        else:
            intfl = '%s_intfile'%(cde)

        # create 4D volume
        cmd = 'fslmerge -t %s'%(os.path.join(outdir,intfl))
        for scn in scans:
            cmd = cmd+' %s'%(scn)

        #if not parallel:
        print 'creating input 4D volume...'
        os.system(cmd)

        #if not parallel:
        print 'loading data'
        data = ni.load(os.path.join(outdir,'%s.nii.gz'%(intfl))).get_data()

    # run voxelwise analysis
    #if not parallel:
    print 'beginning analysis...'
    x,y,z,t_dim = data.shape
    results = np.zeros((x,y,z))
    aff = ni.load(scans[0]).get_affine()

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
        #if not parallel:
        print 'finished %s/%s job clusters'%(xind,x)

    # write image
    outstr = '%s%s'%(os.path.join(outdir,outfl),taskid)
    #if not parallel:
    print 'writing image to %s'%(outstr)
    nimg = ni.Nifti1Image(results,aff)
    ni.save(nimg,outstr)
    outstr = outstr+'.nii'

    # clean up
    data = None
    os.system('rm %s'%(os.path.join(outdir,'%s_*'%(cde))))

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

    scalestr = the string preceding the number indicating atlas resolution in
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
    with a file name indicated by string input


    Outputs a dict where the key is scale and the value is a dataframe
    containing results at that scale
    '''
    save=check_bool(save)
    if save:
        if '_' in save:
            print('WARNING: no _ aloud in save. Removing...')
            nsave=''
            jnk=save.rsplit('_')
            for i in range(len(jnk)):
                nsave=nsave+jnk[i]
            save = nsave

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
            rdf.to_csv(os.path.join(outdir,'%s_scl%s_res%s.csv'%(save,scale,taskid)))

    resdf = pandas.DataFrame(columns =['scale','parcel','measure','value','pvalue'])
    os.system('rm %s'%(os.path.join(outdir,'%s_*'%(cde))))

    return dfz


def id_sig_results(dfz,outdir,outfl='outfl',perc_thr='top',thr_tp='r',thr='',out_tp='samp',res_tp='all',master_ss=False,master_thr='',par=False,parsave=False,taskid=''):
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

    par = If set to the same string as save from the searchlight function, will
    import results spreadsheets generated from searchlight, and will make
    script compatible for command-line based parallelization

    parsave = If true, will also save the spreadsheets imported using par.
    otherwise, spreadsheets will be deleted to conserve space.

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
    par = check_bool(par)
    parsave = check_bool(parsave)

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
    if perc_thr != 'fwe' and perc_thr != 'top':
        if type(thr) != int and type(thr) != float:
            if type(thr) == str:
                try:
                    thr=float(thr)
                except:
                    raise ValueError('thr must be an int or float, unless perc_thr is True')
            else:
                raise ValueError('thr must be an int or float, unless perc_thr is true')
        if thr > 1:
            raise ValueError('invalid value set for thr')
    if perc_thr == 'fwe' and thr_tp != 'p':
        print 'WARNING: because perc_thr set to fwe, change thr_tp to p...'
        thr_tp = 'p'
    if out_tp != 'mast' and outfl=='outfl':
        print 'WARNING: No outfile name specified. Using name %s%s'%(outfl,taskid)

    # wrangle spreadsheets

    if par:
        dfz = parallel_in('isr',outdir,par,parsave,taskid)

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
            master_ss = os.path.join(outdir,'master_ss.csv')
            print 'no master_ss found. Creating new master_ss at %s'%(master_ss)
            mdf = pandas.DataFrame()

    coltest = dfz[dfz.keys()[0]].columns.tolist()
    if res_tp == 'rho' and 'rho' not in coltest:
        raise IOError('res_tp set to rho but no nonparametric results available')
    if res_tp == 'poly' and  len(coltest) < 4:
        raise IOError('res_tp set to poly but no polynomial results available')

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

    res_dict = {'r': 0, 'rho': 2}
    if len(coltest) > 4:
	res_dict.update({'poly': 4})

    if perc_thr:
        comps = 0
        for scl,scf in dfz.iteritems():
            comps = comps + len(scf)
        if res_tp == 'all':
            comps = comps * 3

        if perc_thr == 'fwe':
            thr = (0.05/comps)
            print '%s total comparisons. Setting pvalue threshold to %s'%(comps,thr)
        else:
            vec = []
            for scl,scf in dfz.iteritems():
                for x,y in scf.iterrows():
                    if thr_tp == 'p':
                        if res_tp != 'all':
                            vec.append(y[res_dict[res_tp]+1])
                        else:
                            for k,v in res_dict.iteritems():
                                vec.append(y[v+1])
                    elif thr_tp == 'r' or thr_tp == 'rabs':
                        if res_tp != 'all':
                            vec.append(y[res_dict[res_tp]])
                        else:
                            for k,v in res_dict.iteritems():
                                vec.append(y[v])

            if perc_thr == 'top':
                print 'acquiring top result'
                if thr_tp == 'p':
                    thr = sorted(vec)[1]
                elif thr_tp == 'r' or thr_tp == 'rabs':
                    thr = sorted(vec)[-2]
		print 'threshold < %s'%(thr)
            else:
                frac = int(comps * thr)
                if frac==0:
                    frac=1
                if thr_tp == 'p':
                    thr = sorted(vec)[frac]
                    print 'acquiring top results, p < %s'%(thr)
                else:
                    thr = sorted(vec)[-(frac+1)]
                    print 'acquiring top results, r > %s'%(thr)

    # Extract results
    res_count = 0
    for scl,scf in dfz.iteritems():
        for x,y in scf.iterrows():
            if thr_tp == 'p':
                if res_tp != 'all':
                    if y[(res_dict[res_tp]) + 1] < thr:
                        resdf=update_spreadsheet(resdf,scl,x,y,res_tp)
                        if out_tp != 'samp':
                            if master_thr == thr:
                                mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v='',df_tp='mast',res_count=res_count,taskid=taskid)
                            else:
                                if y[(res_dict[res_tp]) + 1] < master_thr:
                                    mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v='',df_tp='mast',res_count=res_count,taskid=taskid)
                else:
                    for k,v in res_dict.iteritems():
                        if y[v+1] < thr:
                            resdf=update_spreadsheet(resdf,scl,x,y,res_tp,v=v)
                            if out_tp != 'samp':
                                if master_thr == thr:
                                    mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
                                else:
                                    if y[v+1] < master_thr:
                                        mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
            else:
                if res_tp != 'all':
                    if thr_tp == 'rabs':
                        if abs(y[(res_dict[res_tp])]) > thr:
                            resdf=update_spreadsheet(resdf,scl,x,y,res_tp)
                            if out_tp != 'samp':
                                if master_thr == thr:
                                    mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v='',df_tp='mast',res_count=res_count,taskid=taskid)
                                else:
                                    if abs(y[(res_dict[res_tp])]) > thr:
                                        mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v='',df_tp='mast',res_count=res_count,taskid=taskid)
                    else:
                        if y[(res_dict[res_tp])] > thr:
                            resdf=update_spreadsheet(resdf,scl,x,y,res_tp)
                            if out_tp != 'samp':
                                if master_thr == thr:
                                    mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v='',df_tp='mast',res_count=res_count,taskid=taskid)
                                else:
                                    if y[(res_dict[res_tp])] > thr:
                                        mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v='',df_tp='mast',res_count=res_count,taskid=taskid)
                else:
                    for k,v in res_dict.iteritems():
                        if thr_tp == 'rabs':
                            if abs(y[v]) > thr:
                                resdf=update_spreadsheet(resdf,scl,x,y,res_tp,v=v)
                                if out_tp != 'samp':
                                    if master_thr == thr:
                                        mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
                                    else:
                                        if abs(y[v]) > thr:
                                            mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
                        else:
                            if y[v] > thr:
                                resdf = update_spreadsheet(resdf,scl,x,y,res_tp,v=v)
                                print resdf
				if out_tp != 'samp':
                                    if master_thr == thr:
                                        mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)
                                    else:
                                        if y[v] > thr:
                                            mdf,res_count=update_spreadsheet(mdf,scl,x,y,res_tp,v=v,df_tp='mast',res_count=res_count,taskid=taskid)

    # save results

    if out_tp != 'mast':
        outstr = os.path.join(outdir,'%s%s.csv'%(outfl,taskid))
        print 'writing results to spreadsheet at %s'%(outstr)
        resdf.to_csv(outstr)

    if out_tp != 'samp':
        print 'writing results to master spreadsheet at %s'%(master_ss)
        mdf.to_csv(master_ss)

def update_spreadsheet(indf,scl,x,y,res_tp,v='',df_tp='samp',res_count=0,taskid=''):
    '''update spreadsheet with values indexed specifically from searchlight
    generated spreadsheet'''

    res_dict = {'r': 0, 'rho': 2}
    if len(indf.columns.tolist()) > 4:
        res_dict.update({'poly': 4})


    if df_tp == 'samp':
        ind = '%s: %s'%(scl,x)
        cnt = 0
        if ind in indf.index.tolist():
            nind = ind
            while nind in indf.index.tolist():
                cnt = cnt+1
                nind = '%s_%s'%(ind,cnt)
            ind = nind
	# to make compatible with guillimin python
	new_ind = indf.index.tolist()
	new_ind.append(ind)
	indf = indf.reindex(new_ind)

    if df_tp == 'samp':
        if res_tp == 'all':
            indf.ix[ind, 'value'] = y[v]
            indf.ix[ind, 'pvalue'] = y[v+1]
            indf.ix[ind, 'measure'] = res_dict.keys()[res_dict.values().index(v)]
        else:
            indf.ix[ind, 'value'] = y[res_dict[res_tp]]
            indf.ix[ind, 'pvalue'] = y[res_dict[res_tp]+1]
            indf.ix[ind, 'measure'] = res_tp
        indf.ix[ind, 'scale'] = scl
        indf.ix[ind, 'parcel'] = x
	
	return indf

    elif df_tp == 'mast':
	#to make compatible with guillimin python
	new_ind = indf.index.tolist()
	new_ind.append('samp_%s'%(taskid))
	indf = indf.reindex(new_ind)
        cols = ['scale','parcel','value','pvalue','measure']
	for col in cols:
	    ncol = 'res%s_%s'%(res_count,col)
	    ncols = indf.columns.tolist()
	    ncols.append(ncol)
	    indf = indf.reindex(columns=ncols)
            	

        indf.ix['samp_%s'%(taskid),'res%s_scale'%(res_count)] = scl
        indf.ix['samp_%s'%(taskid),'res%s_parcel'%(res_count)] = x
        if res_tp == 'all':
            indf.ix['samp_%s'%(taskid),'res%s_value'%(res_count)] = y[v]
            indf.ix['samp_%s'%(taskid),'res%s_pvalue'%(res_count)] = y[v+1]
            indf.ix['samp_%s'%(taskid),'res%s_measure'%(res_count)] = res_dict.keys()[res_dict.values().index(v)]
        else:
            indf.ix['samp_%s'%(taskid),'res%s_value'%(res_count)] = y[res_dict[res_tp]]
            indf.ix['samp_%s'%(taskid),'res%s_pvalue'%(res_count)] = y[res_dict[res_tp]+1]
            indf.ix['samp_%s'%(taskid),'res%s_measure'%(res_count)] = res_tp
        res_count = res_count + 1

        return indf,res_count

def parallel_in(func,outdir,ipt1,ipt2,ipt3=''):

    if func == 'va':
        idf = pandas.read_csv(ipt1)
        scans = idf[:][idf.columns[0]].tolist()
        pv_vals = idf[:]['pv_vals'].tolist()
	os.system('rm %s'%(ipt1))

        for ind,scan in enumerate(scans):
            if scan[-1] == 'z':
                os.system('cp %s %s/%s_subject%s.nii.gz'%(scan,outdir,ipt2,ind))
            else:
                os.system('cp %s %s/%s_subject%s.nii'%(scan,outdir,ipt2,ind))
        scans = glob(os.path.join(outdir,'%s_*'%(ipt2)))

        return scans,pv_vals

    if func == 'isr':
        dfz = {}
        ssz = glob(os.path.join(outdir,'%s*_res%s.csv'%(ipt1,ipt3)))
        for ss in ssz:
            scl = os.path.split(ss)[1].rsplit('_')[1].rsplit('scl')[1]
            ndf = pandas.read_csv(ss)
            ndf = ndf.drop(ndf.columns[0],axis=1)
            dfz.update({scl: ndf})
            if not ipt2:
                os.remove(ss)

        return dfz

def parallel_out(func,outdir,opt1,opt2,opt3):

    if func == 'dbc':
        odf = pandas.DataFrame(index = opt1,columns = ['pv_vals'])
        for i,sub in enumerate(odf.index.tolist()):
            odf.ix[sub,'pv_vals'] = opt2[i]

    flpth = os.path.join(outdir,'dbc_out%s.csv'%(opt3))
    odf.to_csv(flpth)

    return flpth

def check_bool(var):
    '''for commandline use only. Change string instances of 'False' and 'True'
    into boolean values
    valist = list of variables to convert'''
    
    if var == 'False':
	var = False
    return var

