import os
from glob import glob
import pandas
import random, math
import numpy as np
import nibabel as ni
import scipy.stats.stats as st
import randomstate.prng.mrg32k3a as rnd
import propagation_correlations as wr

def define_bootstrap_sample(ss,subcol,subpth,pv,outpth,resample=False,
                            par=False,rand=False,taskid = '',sample_perc=0.5,
                            num_gps=3):
    """takes an existing sample of subjects from a spreadsheet and 
    creates a subsample, balanced by a variable. Outputs a 
    subsamble membership txt file and a list of paths corresponding 
    to subsample member scans)

    ss = a spreadsheet containing at least a list of subjects and 
    values for a predictor variable

    subcol = string corresponding to the label of the column 
    containing values for the subject IDs. Note these IDs should be 
    consistent with existing file names in subpth. If index, set to 
    'index'

    subpth = a path pointing to the directory containing subject 
    scans

    outpth = desired output directory

    resample =  False = Do not resample data.
                subsamp = Take random subsample of data distributed 
                across pv.*
                * Use sample_perc and num_gps to control parameters
                jackknife = sample using the leave-one-out method
                permute = randomly shuffle sample without 
                replacement
                bootstrap = radnomly shuffle sample with replacement

    pv = string corresponding to the label of the column containing 
    values for the predictor variable

    par = if True, makes compatible for command line based 
    parallelization

    task_id = used to keep track of id in the case of multiple 
    bootstrap samples

    rand = Determine how psuedorandom generator is seeded for the 
    randomization of the sample. Leave as False to use 
    random.shuffle  with a random seed, which will create 
    unreproducible samples and  is not recommended for 
    parallelization. Set as an int to use int as seed for mrg322ka 
    PRNG  (recommended for parallelization or reproducible results)

    sample_perc = For subsamp method, the percent of the total 
    sample to be used in each subsample. Default is 0.5 for 50%. 
    Must be float >0 and <=1

    num_gps = For subsamp method, the number of groups used to 
    balance subsample on predictor variable. Default is 3. Value 
    must be int

    Outputs a list of paths and a vector containing the values for 
    the predictor variable

    ***IMPORTANT NOTE***
    Script has the following assumptions:
    1) scans in subpth directory have ids from subcol in them, but
     the filename of the scan does not /start/ with the id.
    2) There are no more than 1 scan per ID within subpth directory
    """

    par = check_bool(par)
    resample = check_bool(resample)
    rand = check_bool(rand)

    if rand:
        if type(rand) == str:
            rand = int(rand)

    if taskid:
        if type(taskid) == str:
            try:
                taskid=int(taskid)
            except:
                taskid=taskid

    # prep spreadsheet

    if type(num_gps) != int:
        raise TypeError('num_gps must be an integer')

    if sample_perc > 1 or sample_perc <=0:
        raise ValueError('sample_perc must be a float between 0 and 1')

    if ss[-3:] == 'xls' or ss[-4:] == 'xlsx':
        subdf = pandas.ExcelFile(ss).parse('Sheet1')
    elif ss[-3:] == 'csv':
        subdf = pandas.read_csv(ss)
    else:
        raise ValueError(
        'input spreadsheet filetype not recognized. Please use .xls, .xlsx or .csv')

    if subcol != 'index':
        if subcol in subdf.columns.tolist():
            subdf.index = subdf[:][subcol].tolist()
        else:
            raise ValueError('there is no column in the spreadsheet called %s'%
                                                                    (subcol))

    for sub in subdf.index.tolist():
        if len(glob(os.path.join(subpth,'*%s*'%(sub)))) < 1:
            print 'no scan was found for subject %s. Removing from analysis.'%
                                                                        (sub)
            #this was not compatible with guillimin's version of pandas
            #subdf.drop(sub,axis = 0, inplace=True)
            subdf=subdf.drop(sub,axis=0)

    if resample == 'subsamp':
        if pv not in subdf.columns.tolist():
            raise ValueError('there is no column in the spreadsheet called %s'%
                                                                        (pv))
        else:
            subdf = subdf.sort(pv)

    allsubs = subdf.index.tolist()

    # define sample
    if resample == 'subsamp':
        subsamp = subsample(allsubs,num_gps,sample_perc,rand,taskid)
    elif resample == 'jackknife':
        subsamp = jackknife(allsubs,taskid)
    elif resample == 'permute':
        subsamp = permute(allsubs,rand,taskid)
    elif resample == 'bootstrap':
        subsamp = bootstrap(allsubs,rand,taskid)
    elif not resample:
        subsamp = allsubs
    else:
        raise ValueError('value of %s not acceptable for argument resample.'%
                                                                (resample))


    # collect outputs
    scans = []
    for sub in subsamp:
        scn = glob(os.path.join(subpth,'*%s*'%(sub)))[0]
        scans.append(scn)

    if resample == 'subsamp' or resample == 'jackknife':
        pv_vals = []
        for sub in subsamp:
            pv_vals.append(subdf.ix[sub,pv])
    else:
        if pv not in subdf.columns.tolist():
            raise ValueError('there is no column in the spreadsheet called %s'%
                                                                        (pv))
        else:
            pv_vals = subdf[:][pv].tolist()


    # make membership record
    ndf = pandas.DataFrame(index=subsamp)
    if taskid:
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

def randomize_4_montecarlo(seed,jumpseed,jump_mult,rsamp,replace=False):

    intp = type(rsamp[0])
    rs = rnd.RandomState(seed)
    jumpseed = int(jumpseed * jump_mult)
    rs.jump(jumpseed)
    nrsamp = []

    if not replace:
        nrsamp = rs.permutation(rsamp)

    else:
        rand_array = rs.random_integers(0,len(rsamp)-1,len(rsamp))
        for i in rand_array:
            nrsamp.append(rsamp[i])


    return nrsamp

def subsample(allsubs,num_gps,sample_perc,rand,taskid):

   # extract subsample
    if not taskid:
        taskid = 1

    subsamp_n = int(len(allsubs) / num_gps)
    subsamp_dict = {}
    for i in range(1,(num_gps+1)):
        if i == 1:
            subsamp_dict.update({i: allsubs[:subsamp_n]})
        else:
            subsamp_dict.update({i: allsubs[(subsamp_n * (i-1)):(subsamp_n *
            (i))]})

    subsamp = []
    n_per_gp = int(subsamp_n * sample_perc)
    for i,slist in subsamp_dict.iteritems():
        if not rand:
            rlist = random.shuffle(slist)
        else:
            rlist = randomize_4_montecarlo(rand,taskid,(i+1),slist)
        for x in range(n_per_gp):
            subsamp.append(rlist[x]) 

    return subsamp

def jackknife(sample,taskid=None):

    if not taskid:
        taskid=1

    if taskid >= len(sample):
        raise ValueError(
            'With jackknife, only %s samples possible. Aborting...'%
                                                    (len(allsubs)))
    else:
        dropper = sample[taskid - 1]
        sample.remove(dropper)

        return sample

def permute(sample,rand,taskid=None):
    
    if not taskid:
        taskid=1

    if not rand:
        rlist = random.shuffle(sample)
    else:
        rlist = randomize_4_montecarlo(rand,taskid,1,sample)

    return rlist

def bootstrap(sample,rand,taskid=None):

    if not taskid:
        taskid=1

    if not rand:
        rlist = []
        for i in range(len(sample)):
            j = random.choice(sample)
            rlist.join(j)

    else:
        rlist = randomize_4_montecarlo(rand,taskid,1,sample,replace=True)

    return rlist

def convert_r2t(r,samp_df):
    t = math.sqrt((samp_df * r**2) / (1 - r**2))
    if r < 0:
        t = (t*-1)

    return t

def voxelwise_analysis(scans,pv_vals,outfl,outdir,out_tp='r',nonpar=False,
                       taskid='',parallel=False,parin='',indata=False,
                       intermed=False):
    """Given a list of 3D volumes and a corresponding list of 
    values for a predictor variable, run voxelwise correlations.

    scans = a list of paths corresponding to 3D nifti volumes

    pv_vals = a list of values corresponding to subject values, in 
    the same order as scans

    outfl = string to be used to identify the outfile map (NOTE, 
    taskid will be automatically appended to this string)

    outdir = the path to the directory where the outfiles will be 
    written to

    out_tp = the type of map to be generated. 'r' will generate a 
    voxelwise rmap, whereas 't' will generate a voxelwise tmap

    nonpar = Set to True to set voxelwise analysis from pearson to 
    spearman

    indata = In case 4D data is already available in an existing 
    variable, data can be specified here. Or, if an intermed file 
    exists, simply add the path here.

    taskid = used to keep track of id in the case of multiple 
    bootstrap samples

    parallel = if True, script will copy scans into working 
    directory to allow for multiple concurrent processes. Will also 
    make script compatible for command line based parallelization.

    parin = input path for parallelization

    intermed = if True, script will not delete the 4D volume used 
    to run the voxelwise analysis


    Outputs a path pointing to the newly created tmap


    WARNING: As of now, this script is extremely computationally 
    intensive for a local machine. Running a subsample of 135 
    subjects required 5 GB memory, used 100% of my CPU at times, 
    and took well over 10 minutes...

    NOTE: As of now, script does not regress out confounding 
    variables or mask analysis. Also, much more efficient voxelwise 
    analysis code is on the way...
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
            if not parallel:
                print 'loading data...'
            data=ni.load(indata).get_data()
    else:
        if intermed:
            intfl = 'intfile%s'%(taskid)
        else:
            intfl = '%s_intfile'%(cde)

        # create 4D volume
#        cmd = 'fslmerge -t %s'%(os.path.join(outdir,intfl))
#        for scn in scans:
#            cmd = cmd+' %s'%(scn)

#        print 'creating input 4D volume...'
#       os.system(cmd)

#        print 'loading data'
#        data = ni.load(os.path.join(outdir,'%s.nii.gz'%(intfl))).get_data()

        print 'create 4D volume'

        to_make = []
        for scn in scans:
            nscn = ni.load(scn).get_data()
            to_make.append(nscn)
        data = np.concatenate([aux[..., None] for aux in to_make], axis=3)



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
    os.system('rm %s'%(os.path.join(outdir,'%s_*'%(cde))))

    return outstr

def spatial_correlation_searchlight_from_NIAK_GLMs(indir,cov_img,outdir,
                                                   contrast,templ_str,
                                                   scalestr='scale',
                                                   norm=False,xfm=False,
                                                   eff='t',poly=1,taskid='',
                                                   save=False):
    '''
    Given a directory containing a) NIAK glm.mat files at multiple 
    scales, b) atlases at the same resolutions, and a covariate 
    image, this function will do the following:
        1) Convert covariate image to resolution of each atlas
        2) Perform spatial correlations between a) average 
        connectivity matrix seed from every region at every scale, 
        b) the covariance image
        3) Output a dictionary summarizing the results, where key 
        = resolution and value = a dataframe containing statistics 
        for each test


    indir = path to directory containing NIAK glm.mat files and 
    atlases.

    cov_img = path to covariate image

    outdir = path to desired output directory

    contrast = label of the desired contrast to extract average 
    connectivity values from glm.mat file

    templ_str = search string to locate atlases

    scalestr = the string preceding the number indicating atlas 
    resolution in the glm and atlas paths

    norm = If set to a path, will normalize cov_img to target image 
    that path points to. Uses flirt with nearest neighbour 
    interpolation. If False, no normalization will occur. 
    Acceptable inputs are nii or xfm files

    eff = if set to 'r', will take values within the 'eff' 
    structure of the glm.mat file. If set to 't' will take values 
    within the 'ttest' structure of the glm.mat file

    poly = If int > 1, results will include test statistics modeled 
    with a n-order polynomial, where n = poly

    taskid = used to keep track of id in the case of multiple 
    bootstrap samples

    save = If set to a string, will write results from each 
    resolution to a spreadsheet with a file name indicated by 
    string input


    Outputs a dict where the key is scale and the value is a 
    dataframe containing results at that scale
    '''
    save = check_bool(save)
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
	    # get rid of nans if necessary
	    data = ni.load(cov_img).get_data()
	    if not pandas.notnull(data[0][0][0]):
	        nimg = os.path.join(outdir,'%s_ncov'%(cde))
	        os.system('fslmaths %s -nan %s'%(cov_img,nimg))
	        cov_img = '%s.nii.gz'%(nimg)	

        if xfm:
	        os.system(
                'flirt -interp nearestneighbour -in %s -ref %s -applyxfm -init %s -out %s/%s_rtmap'%
                                                (cov_img,norm,xfm,outdir,cde))
      	    else:
                os.system(
                    'flirt -interp nearestneighbour -in %s -ref %s -out %s/%s_rtmap'%
                                                    (cov_img,norm,outdir,cde))
        cov_img = os.path.join(outdir,'%s_rtmap.nii.gz'%(cde))

    glmz = glob(os.path.join(indir,'glm*.mat'))
    dfz = {}
    for glm in glmz:
    ####################the 2 lines below suck######################
        scale = int(os.path.split(glm)[1].rsplit('_%s'%
                                                (scalestr))[1].rsplit('.')[0])
        scale_templ = glob(os.path.join(indir,'%s*%s.*'%(templ_str,scale)))[0]
        df, rdf, scalar = wr.get_parcelwise_correlation_map_from_glm_file(
                outdir,glm,scale_templ,scale,contrast,eff=eff,cov_msk='',
                cov_img=cov_img,conndf = '',poly=poly)
        dfz.update({scale: rdf})
        if save:
            rdf.to_excel(os.path.join(outdir,'%s_scl%s_res%s.xls'
                                                        %(save,scale,taskid)))

    resdf = pandas.DataFrame(columns =[
                                'scale','parcel','measure','value','pvalue'])
    os.system('rm %s'%(os.path.join(outdir,'%s_*'%(cde))))

    return dfz


def id_sig_results(dfz,outdir,outfl='outfl',perc_thr='top',thr_tp='r',thr='',
                   out_tp='samp',res_tp='all',master_ss=False,master_thr='',
                   par=False,parsave=False,taskid=''):
    '''given a dict outputted from the searchlight funcion, and 
    thresholding information, this function will save the top 
    results from the searchlight into a run specific spreadsheet 
    and/or a master spreadsheet (across bootstrap samples)

    dfz = a dict outputted from the searchlight function where the 
    key is atlas resolution and value is a specific dataframe of 
    results

    outdir = the desired output directory

    outfl = the name of the  output file (NOTE: taskid will be
    automatically appended to this)

    perc_thr: 'perc' = set threshold by top percentage of values 
              (values determined by thr_tp)
              'top' = keep only the top result (value determined 
              by res_tp)
              'fwe' = threshold based on # of comparisons 
              (bonferonni)

    out_tp =  'samp' = only create output spreadsheet for the input
              'mast' = only write results to a master spreadsheet
              'both' = do both

    thr_tp =  'p' = set threshold by pvalue
              'r' = set threshold by test statistic
              'rabs' = set threshold by absolute value of test 
               statistic

    thr = desired significance threshold if thr_tp set to 'p','r', 
    or if perc_thr set to 'perc'

    res_tp =  'r' = only report values from standard pearson r 
               results
              'rho' = only report values from spearman r results
              'poly' = only report values from polynomial results
              'all' = report values from all three types of tests

    master_ss = path to the master spreadsheet. If no path is 
    given, a master spreadsheet will be created. Will be ignored 
    if out_tp is set to 'samp'

    master_thr = threshold for master spreadsheet results. If 
    blank, will use value from thr. Will be ignored if out_tp set 
    to 'samp'. NOTE: master_thr MUST be equal to or more 
    conservative than thr

    par = If set to the same string as save from the searchlight 
    function, will import results spreadsheets generated from 
    searchlight, and will make script compatible for command-line 
    based parallelization

    parsave = If True, will also save the spreadsheets imported 
    using par. otherwise, spreadsheets will be deleted to conserve 
    space.

    taskid = used to keep track of id in the case of multiple 
    bootstrap samples

    WARNING: To avoid giant files and potentially breaking things, 
    consider setting thr to very conservative stats, especially if 
    res_tp = 'all'. Even 'fwe' will often yield >5% of results 
    significant. Using perc and low values, or r and a high thr, 
    is recommended.

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
                raise ValueError('thr must be an int or float, unless perc_thr is True')
        if thr > 1:
            raise ValueError('invalid value set for thr')
    if perc_thr == 'fwe' and thr_tp != 'p':
        print 'WARNING: because perc_thr set to fwe, changing thr_tp to p...'
        thr_tp = 'p'
    if out_tp != 'mast' and outfl=='outfl':
        print 'WARNING: No outfile name specified. Using name %s%s'%
                                                    (outfl,taskid)

    # wrangle spreadsheets

    if par:
        dfz = parallel_in('isr',outdir,par,parsave,taskid)

    resdf = pandas.DataFrame(columns =[
            'scale','parcel','measure','value','pvalue'])

    if out_tp != 'samp':
        if os.path.isfile(master_ss):
            if master_ss[-3:] == 'xls' or master_ss[-4:] == 'xlsx':
                mdf = pandas.ExcelFile(master_ss).parse('Sheet1')
            elif master_ss[-3:] == 'csv':
                mdf = pandas.read_csv(master_ss)
            else:
                raise ValueError('input spreadsheet filetype not recognized')
        else:
            master_ss = os.path.join(outdir,'master_ss.xls')
            print 'no master_ss found. Creating new master_ss at %s'%
                                                                (master_ss)
            mdf = pandas.DataFrame()

    coltest = dfz[dfz.keys()[0]].columns.tolist()
    if res_tp == 'rho' and 'rho' not in coltest:
        raise ValueError(
                'res_tp set to rho but no nonparametric results available')
    if res_tp == 'poly' and len(coltest) < 8:
        raise ValueError(
                'res_tp set to poly but no polynomial results available')

    # determine thresholds    
    if out_tp != 'samp':
        if master_thr == '':
            master_thr = thr

        if thr_tp == 'r' or thr_tp == 'rabs':
            if master_thr < thr:
                raise ValueError(
                        'master_thr must be more conservative than thr')
        else:
            if master_thr > thr:
                raise ValueError(
                        'master_thr must be more conservative than thr')

    res_dict = {'r': 0, 'rho': 2}
    if len(coltest) > 8:
        res_dict.update({'poly': 4})

    if perc_thr:
        comps = 0
        for scl,scf in dfz.iteritems():
            comps = comps + len(scf)
        if res_tp == 'all':
            comps = comps * 3

        if perc_thr == 'fwe':
            thr = (0.05/comps)
            print '%s total comparisons. Setting pvalue threshold to %s'%
                                                                    (comps,thr)
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
            else:
                frac = int(comps * thr)
                if frac==0:
                    frac=1
                if thr_tp == 'p':
                    thr = sorted(vec)[frac]
                    print 'acquiring top results, p < %s'%(thr)
                else:
                    if thr == 1.0:
                        thr = sorted(vec)[-(frac)]
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
                        update_spreadsheet(resdf,scl,x,y,res_tp)
                        if out_tp != 'samp':
                            if master_thr == thr:
                                res_count=update_spreadsheet(
                                        mdf,scl,x,y,res_tp,v='',df_tp='mast',
                                        res_count=res_count,taskid=taskid)
                            else:
                                if y[(res_dict[res_tp]) + 1] < master_thr:
                                    res_count=update_spreadsheet(
                                            mdf,scl,x,y,res_tp,v='',df_tp='mast',
                                            res_count=res_count,taskid=taskid)
                else:
                    for k,v in res_dict.iteritems():
                        if y[v+1] < thr:
                            update_spreadsheet(resdf,scl,x,y,res_tp,v=v)
                            if out_tp != 'samp':
                                if master_thr == thr:
                                    res_count=update_spreadsheet(
                                            mdf,scl,x,y,res_tp,v=v,
                                            df_tp='mast',res_count=res_count,
                                            taskid=taskid)
                                else:
                                    if y[v+1] < master_thr:
                                        res_count=update_spreadsheet(
                                                mdf,scl,x,y,res_tp,v=v,
                                                df_tp='mast',
                                                res_count=res_count,
                                                taskid=taskid)
            else:
                if res_tp != 'all':
                    if thr_tp == 'rabs':
                        if abs(y[(res_dict[res_tp])]) > thr:
                            update_spreadsheet(resdf,scl,x,y,res_tp)
                            if out_tp != 'samp':
                                if master_thr == thr:
                                    res_count=update_spreadsheet(
                                            mdf,scl,x,y,res_tp,v='',
                                            df_tp='mast',res_count=res_count,
                                            taskid=taskid)
                                else:
                                    if abs(y[(res_dict[res_tp])]) > thr:
                                        res_count=update_spreadsheet(
                                                mdf,scl,x,y,res_tp,v='',
                                                df_tp='mast',
                                                res_count=res_count,
                                                taskid=taskid)
                    else:
                        if y[(res_dict[res_tp])] > thr:
                            update_spreadsheet(resdf,scl,x,y,res_tp)
                            if out_tp != 'samp':
                                if master_thr == thr:
                                    res_count=update_spreadsheet(
                                            mdf,scl,x,y,res_tp,v='',
                                            df_tp='mast',res_count=res_count,
                                            taskid=taskid)
                                else:
                                    if y[(res_dict[res_tp])] > thr:
                                        res_count=update_spreadsheet(
                                                mdf,scl,x,y,res_tp,v='',
                                                df_tp='mast',
                                                res_count=res_count,
                                                taskid=taskid)
                else:
                    for k,v in res_dict.iteritems():
                        if thr_tp == 'rabs':
                            if abs(y[v]) > thr:
                                update_spreadsheet(resdf,scl,x,y,res_tp,v=v)
                                if out_tp != 'samp':
                                    if master_thr == thr:
                                        res_count=update_spreadsheet(
                                                mdf,scl,x,y,res_tp,v=v,
                                                df_tp='mast',
                                                res_count=res_count,
                                                taskid=taskid)
                                    else:
                                        if abs(y[v]) > thr:
                                            res_count=update_spreadsheet(
                                                    mdf,scl,x,y,res_tp,v=v,
                                                    df_tp='mast',
                                                    res_count=res_count,
                                                    taskid=taskid)
                        else:
                            if y[v] > thr:
                                update_spreadsheet(resdf,scl,x,y,res_tp,v=v)
                                if out_tp != 'samp':
                                    if master_thr == thr:
                                        res_count=update_spreadsheet(
                                                mdf,scl,x,y,res_tp,v=v,
                                                df_tp='mast',
                                                res_count=res_count,
                                                taskid=taskid)
                                    else:
                                        if y[v] > thr:
                                            res_count=update_spreadsheet(
                                                    mdf,scl,x,y,res_tp,v=v,
                                                    df_tp='mast',
                                                    res_count=res_count,
                                                    taskid=taskid)

    # save results

    if out_tp != 'mast':
        outstr = os.path.join(outdir,'%s_thrres%s.xls'%(outfl,taskid))
        print 'writing results to spreadsheet at %s'%(outstr)
        resdf.to_excel(outstr)

    if out_tp != 'samp':
        print 'writing results to master spreadsheet at %s'%(master_ss)
        mdf.to_excel(master_ss)

def update_spreadsheet(indf,scl,x,y,res_tp,
                       v='',df_tp='samp',
                       res_count=0,taskid=''):
    '''update spreadsheet with values indexed specifically from 
    searchlight generated spreadsheet'''

    res_dict = {'r': 0, 'rho': 2}
    if res_tp == 'all' or res_tp == 'poly':
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

    if df_tp == 'samp':
        if res_tp == 'all':
            indf.ix[ind, 'value'] = y[v]
            indf.ix[ind, 'pvalue'] = y[v+1]
            indf.ix[ind, 'measure'] = res_dict.keys()[
                                                res_dict.values().index(v)]
        else:
            indf.ix[ind, 'value'] = y[res_dict[res_tp]]
            indf.ix[ind, 'pvalue'] = y[res_dict[res_tp]+1]
            indf.ix[ind, 'measure'] = res_tp
        indf.ix[ind, 'scale'] = scl
        indf.ix[ind, 'parcel'] = x

    elif df_tp == 'mast':
        indf.ix['samp_%s'%(taskid),'res%s_scale'%(res_count)] = scl
        indf.ix['samp_%s'%(taskid),'res%s_parcel'%(res_count)] = x
        if res_tp == 'all':
            indf.ix['samp_%s'%(taskid),'res%s_value'%(res_count)] = y[v]
            indf.ix['samp_%s'%(taskid),'res%s_pvalue'%(res_count)] = y[v+1]
            indf.ix[
                    'samp_%s'%(taskid),
                    'res%s_measure'%(res_count)] = res_dict.keys()[
                                                res_dict.values().index(v)]
        else:
            indf.ix[
                    'samp_%s'%(taskid),
                    'res%s_value'%(res_count)] = y[res_dict[res_tp]]
            indf.ix[
                    'samp_%s'%(taskid),
                    'res%s_pvalue'%(res_count)] = y[res_dict[res_tp]+1]
            indf.ix[
                    'samp_%s'%(taskid),
                    'res%s_measure'%(res_count)] = res_tp
        res_count = res_count + 1

        return res_count

def parallel_in(func,outdir,ipt1,ipt2,ipt3=''):
    '''handles inputs for the various functions when in the case 
    of parallelization
    '''

    if func == 'va':
        idf = pandas.ExcelFile(ipt1).parse('Sheet1')
        scans = idf.index.tolist()
        pv_vals = idf[:]['pv_vals'].tolist()
        #os.system('rm %s'%(ipt1))

        for ind,scan in enumerate(scans):
            if len(str(ind)) == 1:
		ind = '00%s'%(ind)
	    elif len(str(ind)) == 2:
		ind = '0%s'%(ind)
	    if scan[-1] == 'z':
                os.system('cp %s %s/%s_subject%s.nii.gz'%
                                            (scan,outdir,ipt2,ind))
            else:
                os.system('cp %s %s/%s_subject%s.nii'%(scan,outdir,ipt2,ind))
        scans = sorted(glob(os.path.join(outdir,'%s_*'%(ipt2))))

        return scans,pv_vals

    if func == 'isr':
        dfz = {}
        ssz = glob(os.path.join(outdir,'%s*_res%s.xls'%(ipt1,ipt3)))
        for ss in ssz:
            scl = os.path.split(ss)[1].rsplit('_')[1].rsplit('scl')[1]
            ndf = pandas.ExcelFile(ss).parse('Sheet1')
            dfz.update({scl: ndf})
            if not ipt2:
                os.remove(ss)

        return dfz

    if func == 'coi':
        if ipt1[-1] == 'v':
            idf = pandas.read_csv(ipt1)
            idf.index = idf[:][idf.columns[0]].tolist()
        else:
            idf = pandas.ExcelFile(ipt1).parse('Sheet1')

        return idf


def parallel_out(func,outdir,opt1,opt2,opt3):
    '''handles outputs for the various functions when in the case 
    of parallelization
    '''

    if func == 'dbc':
        odf = pandas.DataFrame(index = opt1,columns = ['pv_vals'])
        for i,sub in enumerate(odf.index.tolist()):
            odf.ix[sub,'pv_vals'] = opt2[i]

    flpth = os.path.join(outdir,'dbc_out%s.xls'%(opt3))
    odf.to_excel(flpth)

    return flpth

def check_bool(var):
    '''for commandline use only. Change string instances of 'False' and 'True'
    into boolean values
    '''

    if var == 'False':
        var = False
    return var

def collect_results(ss_dir,ss_str,ss_ext,outdir,outfl='',thr_tp='r',
                    resample=False,permtest='',summary = False):
    '''Given result spreadsheets generated by id_sig_results, will
    summarize data from all spreadsheets into single summary 
    spreadsheets. Function will automatically handle results using 
    multiple statistics, and will output separate results files 
    for them.

    ss_dir = directory where results spreadsheets are stored

    ss_str = search string unique to spreadsheets. Assumes that 
    search string is the beginning of the filenames and that, with 
    the extension, ssdir/ss_str*.ss_ext will find only the results 
    spreadsheets

    ss_ext = String extention for results files. Will also 
    determine extension of outfiles. No '.' needed. Example: 'xls'

    outdir = Directory where you wish the results spreadsheets to 
    be outputted

    outfl = string to label outfiles

    thr_tp =  'r' = rank results using coefficient (r, rho, etc)
              'rabs' = rank results using absolute value of 
              coefficient
              'p' = ranks results using pvalue

    resample = if not False, will assess distribution of test 
    statistics or pvalues (set with thr_tp) and output the lower 
    percentage of this distribution, the value entered into 
    resample (between 0 and 1, should be something like 0.05 for 
    a alpha of 0.05) divided by two (i.e. if resample is 0.05, you 
    would get the 2.5% and 97.5% values).

    permtest = If not False, use random sampling results to assess 
    where in sampling distribution a value falls, and whether it 
    is below an alpha threshold. Input can be a pandas dataframe 
    or a path to a spreadsheet that has statistics from one sample 
    (probably the original sample). Will use input of resample as 
    the critical value.

    summary = If True, will output an excel measure that will 
    contain information about top results across all bootstrap 
    samples. This is somewhat time consuming and the results are 
    not very meaningful.


    Results spreadsheet will contain information on the following:
    top_hits =     number of times a seedmap had the strongest 
                   result of its batch
    appearances =  number of times a seedmap appeared in the 
                   thresholded results
    stat_mean =    the average coefficient for this seedmap across 
                   batches
    stat_sd =      the SD of the coefficient for this seedmap 
                   across boatches
    tophits_sign = the direction of the effect, either pos, neg, 
                   or both

    If certain flags are thrown, may also contain:
    tophit_sign =     indicates whether tophit has a positive or 
                      negative sign
    upper_lim_stat =  the upper limit of the resample stat 
                      distribution, as set by resample
    lower_lim_stat =  the lower limit of the resample stat 
                      distribution, as set by resample
    resample_pvalue = the probability the statistic is not by 
                      chance
    resample_sig =    whether the statistic is significantly not 
                      likely to be chance    

    '''
    resample = check_bool(resample)
    summary = check_bool(summary)

    if resample:
        if type(resample) == str:
            resample = float(resample)
        if resample > 1 or resample < 0:
            raise ValueError('resample must either be a value between 1 and 0, or must be set to False')

    if type(permtest) == pandas.core.frame.DataFrame:
        compdf = permtest
        permtest = True
    else:
        if permtest:
            if type(permtest) == str:
                if permtest[-1] == 'v':
                    compdf = pandas.read_csv(permtest)
                    compdf.index = compdf[:][compdf.columns[0]].tolist()
                else:
                    print permtest
                    compdf = pandas.ExcelFile(permtest).parse('Sheet1')
            else:
                raise ValueError('input for permtest invalid. Please pass a pandas dataframe or a path to a spreadsheet')

    compdf = remove_redundancy_labels(compdf)

    if thr_tp == 'r' or 'rabs':
        vcol = 'value'
    elif thr_tp == 'p':
        vcol = 'pvalue'

    ssz = sorted(glob(os.path.join(ss_dir,'%s_thrres*.%s'%(ss_str,ss_ext))))

    # concat all frames for information and top results
    print 'concatenating frames'
    framez = []
    for ss in ssz:
        if ss_ext[0] == 'x':
            df = pandas.ExcelFile(ss).parse('Sheet1')
        else:
            df = pandas.read_csv(ss)
            df.index = df[:][df.columns[0]].tolist()
        framez.append(df)
    bigdf = pandas.concat(framez)
    st_typez = pandas.unique(bigdf[:]['measure'])
    sumcolz = []
    coltpz = ['_top','_topreg']
    if thr_tp == 'rabs':
        coltpz.append('_sign')
    for tp in st_typez:
        for ctp in coltpz:
            col = '%s%s'%(tp,ctp)
            sumcolz.append(col)
    sumdf = pandas.DataFrame(np.zeros((1,len(sumcolz))),columns = sumcolz)
    bigdf = bigdf.sort('measure')

    # address redundancy labeling from id_sig_results
    print 'removing redundancy labeling...'
    bigdf = remove_redundancy_labels(bigdf)

    # get useless summary measures
    if summary:
        print 'creating summary measures'
        for tp in st_typez:
            for i,indz in enumerate(bigdf.iterrows()):
                if indz[1]['measure'] == tp:
                    if i == 0:
                        if thr_tp == 'rabs':
                            sumdf.ix[
                                    sumdf.index[0],
                                    '%s_top'%(tp)] = abs(indz[1][vcol])
                        else:    
                            sumdf.ix[
                                    sumdf.index[0],
                                    '%s_top'%(tp)] = indz[1][vcol]
                        sumdf.ix[
                                sumdf.index[0],
                                '%s_topreg'%(tp)] = indz[0]
                    else:
                        if thr_tp == 'r':
                            if indz[1][vcol] > sumdf.ix[
                                                        sumdf.index[0],
                                                        '%s_top'%(tp)]:
                                sumdf.ix[
                                        sumdf.index[0],
                                        '%s_top'%(tp)] = indz[1][vcol]
                                sumdf.ix[
                                        sumdf.index[0],
                                        '%s_topreg'%(tp)] = indz[0]  
                        elif thr_tp == 'rabs':
                            if abs(indz[1][vcol]) > sumdf.ix[
                                                            sumdf.index[0],
                                                            '%s_top'%(tp)]:
                                sumdf.ix[
                                        sumdf.index[0],
                                        '%s_top'%(tp)] = abs(indz[1][vcol])
                                sumdf.ix[
                                        sumdf.index[0],
                                        '%s_topreg'%(tp)] = indz[0]
                                if indz[1][vcol] < 0:
                                    sumdf.ix[
                                            sumdf.index[0],
                                            '%s_sign'%(tp)] = 'neg'
                                else:
                                    sumdf.ix[
                                            sumdf.index[0],
                                            '%s_sign'%(tp)] = 'pos'
                        else:
                            if indz[1][vcol] < sumdf.ix[
                                                        sumdf.index[0],
                                                        '%s_top'%(tp)]:
                                 sumdf.ix[
                                        sumdf.index[0],
                                        '%s_top'%(tp)] = indz[1][vcol]
                                 sumdf.ix[
                                        sumdf.index[0],
                                        '%s_topreg'%(tp)] = indz[0]
        if ss_ext[0] == 'x':
            sumdf.to_excel(os.path.join(outdir,'summary.xls'))
        else:
            sumdf.to_csv(os.path.join(outdir,'summary.csv'))

    # create main result dataframe
    res_regz = pandas.unique(bigdf.index.tolist())

    if len(st_typez) == 1:
        cols = ['top_hits','appearances']
        cols.append('mean_%s'%(st_typez[0]))
        cols.append('sd_%s'%(st_typez[0]))
        if thr_tp == 'rabs':
            cols.append('tophit_sign')
        if resample:
            mcols = ['upper_lim_%s'%(st_typez[0]),'lower_lim_%s'%
                                                            (st_typez[0])]
            for mc in mcols:
                cols.append(mc)
            if permtest:
                cols.append('resample_pvalue')
                cols.append('resample_sig')
                cols.append('test_value')
        rdf = pandas.DataFrame(index = res_regz, columns = cols)
        
        # get appearances and confidence statistics
        print 'getting appearances and confidence statistics'
        for lab in res_regz:
            ind = bigdf.ix[lab,vcol]
            synthesize_results(res_regz, bigdf, rdf, vcol,
                                st_typez[0],resample,compdf)

    else:
        rdfz = {}
        for tp in st_typez:
            rind = []
            for lab in res_regz:
                ind = bigdf.ix[lab,'measure']
                if type(ind) == pandas.core.series.Series:
                    if tp in ind.tolist():
                        rind.append(lab)
                else:
                    if ind == tp:
                        rind.append(lab)
            cols = ['top_hits','appearances','mean_%s'%(tp),'sd_%s'%(tp)]
            if thr_tp == 'rabs':
                cols.append('tophit_sign')
            if resample:
                mcols = ['upper_lim_%s'%(tp),'lower_lim_%s'%(tp)]
                for mc in mcols:
                    cols.append(mc)
                if permtest:
                    cols.append('resample_pvalue')
                    cols.append('resample_sig')
                    cols.append('test_value')
            rdf = pandas.DataFrame(index = rind,columns = cols)
            rdfz.update({tp: rdf})

        # get appearances and confidence statistics
        print 'getting appearances and confidence statistics'

        for tp,rdf in rdfz.iteritems():
            if type(compdf) == pandas.core.frame.DataFrame:
                tstdf = compdf.loc[compdf['measure'] == tp]
            ndf = bigdf.loc[bigdf['measure'] == tp]
            synthesize_results(rdf.index.tolist(), ndf, rdf, 
                                vcol, tp, resample, tstdf)

    # extract top hits
    bigdf = None #don't need it any more and takes up a lot of memory
    print 'getting top hits...'
    for ss in ssz:
        if ss_ext[0] == 'x':
            df = pandas.ExcelFile(ss).parse('Sheet1')
        else:
            df = pandas.read_csv(ss)
            df.index = df[:][df.columns[0]].tolist()
        if len(st_typez) == 1: 
            extract_top_hit(rdf,df,vcol,thr_tp)
        else:
            df = remove_redundancy_labels(df)
            for tp,rdf in rdfz.iteritems():
                ndf = df.loc[df['measure'] == tp]
                if len(ndf) > 0:
                    extract_top_hit(rdf,ndf,vcol,thr_tp)

    # save results
    
    if len(st_typez) == 1:
        if ss_ext[0] == 'x':
            rdf.to_excel(os.path.join(outdir,'%s_finalres.xls'%(outfl)))
            print 'spreadsheet being written to %s_finalres.xls'%(outfl)
        else:
            rdf.to_csv(os.path.join(outdir,'%s_finalres.csv'%(outfl)))
            print 'spreadsheet being written to %s_finalres.csv'%(outfl)
    else:
        for tp,rdf in rdfz.iteritems():
            if ss_ext[0] == 'x':
                rdf.to_excel(os.path.join(outdir,'%s_finalres%s.xls'%
                                                                (outfl,tp)))
                print 'spreadsheet being written to %s_finalres%s.xls'%
                                                                (outfl,tp)
            else:
                rdf.to_csv(os.path.join(outdir,'%s_finalres%s.csv'%
                                                                (outfl,tp)))  
                print 'spreadsheet being written to %s_finalres%s.csv'%
                                                                (outfl,tp)          

def remove_redundancy_labels(df):
    '''for its own reasons, id_sig_results creates separate labels for
    seedmaps that have "significant" results from multiple types of tests. This
    function gets rid of the redundancy labels so collect_results can run
    properly

    df = input dataframe

    Outputs a dataframe with redudancy labels removed from the index'''

    nindx = []
    for ind in df.index.tolist():
        if '_' in ind:
            nind,_ = ind.split('_')
            nindx.append(nind)
        else:
            nindx.append(ind)
    df.index = nindx 

    return df

def synthesize_results(res_regz, bigdf, rdf, vcol,
                        st_type,resample,tstdf=False):
    for lab in res_regz:
        ind = bigdf.ix[lab,vcol]
        if type(ind) == pandas.core.series.Series:
            rdf.ix[lab,'appearances'] = len(ind)
            rdf.ix[lab,'mean_%s'%(st_type)] = ind.mean()
            rdf.ix[lab,'sd_%s'%(st_type)] = ind.std()
            if resample:
                p1 = int((resample/2)*len(ind))
                p2 = int((1-(resample/2))*len(ind))
                val1 = sorted(ind)[p1]
                val2 = sorted(ind)[p2]
                rdf.ix[lab,'upper_lim_%s'%(st_type)] = val2
                rdf.ix[lab,'lower_lim_%s'%(st_type)] = val1
                if type(tstdf) == pandas.core.frame.DataFrame:
                    p = resample_test(ind,lab,tstdf,vcol)
                    if p != []:
                        rdf.ix[lab,'resample_pvalue'] = p
                        rdf.ix[lab,'test_value'] = tstdf.ix[lab,vcol]
                        if p < resample:
                            rdf.ix[lab,'resample_sig'] = 1
                        else:
                            rdf.ix[lab,'resample_sig'] = 0


        else:
            rdf.ix[lab,'appearances'] = 1
            rdf.ix[lab,'mean_%s'%(st_type)] = ind
            rdf.ix[lab,'sd_%s'%(st_type)] = np.nan
            if resample:
                rdf.ix[lab,'upper_lim_%s'%(st_type)] = np.nan
                rdf.ix[lab,'lower_lim_%s'%(st_type)] = np.nan

def resample_test(distrib,lab,tstdf,vcol):
    slicer = []
    if lab in tstdf.index.tolist():
        tstval = tstdf.ix[lab,vcol]
        for val in sorted(distrib):
            if vcol == 'value':
                if val > tstval:
                    slicer.append(val)
            elif vcol == 'pvalue':
                if val < tstval:
                    slicer.append(val)
            else:
                raise IOError('vcol not entered properly')
        
        p = len(slicer)/float(len(distrib))

    else:
        p = []

    return p

def extract_top_hit(rdf,df,vcol,thr_tp):

    df = df.sort(vcol)
    if thr_tp == 'r':
        hit = df.index[-1]
    elif thr_tp == 'p':
        hit = df.index[0]
    elif thr_tp == 'rabs':
        for sub in df.index.tolist():
            df.ix[sub,'abs'] = df.ix[sub,'value']
        df = df.sort('abs')
        hit = df.index[-1]
        if df.ix[hit,vcol] > 0:
            hitsign = 'pos'
        else:
            hitsign = 'neg'
    
    prev = rdf.ix[hit,'top_hits']
    if str(prev) == 'nan':
        rdf.ix[hit,'top_hits'] = 1
        if thr_tp == 'rabs':
            rdf.ix[hit,'tophit_sign'] = hitsign
    else:
        rdf.ix[hit,'top_hits'] = (prev+1)
        if thr_tp == 'rabs':
            if rdf.ix[hit,'tophit_sign'] != hitsign:
                rdf.ix[hit,'tophit_sign'] = 'both'

#def 

def create_output_images(resdf,scldir,sclstr,outdir,outstr,input_tp,
                        output_tp='resamp',par=False):
    '''Using results from function collect_results as input, will 
    generate parcel-wise brain images as a spatial representation 
    of results.

    resdf = a pandas dataframe generated from function 
    collect_results. Will be ignored if par is not False.

    scldir = a path pointing to the directory containing atlases 
    for each resolution represented in the results.

    sclstr = string unique to the atlas such that 
    [scldir]/[scstr][scale].* will successfully find each atlas

    outdir = path to desired output directory

    outstr = label to ID output images

    input_tp = as of now, script will only create outputs for one 
    "measure" at a time. Indicate the type of measure (i.e. a 
    valid value within the measure column of the results). 
    Example: 'rho'

    output_tp = resamp = Will only output resampled p-value image
                cross = Will output appearances and top hits, and 
                average stat
                ci = Will only output average stat
                all = all output images will be generated

    par = Path pointing to an input results spreadsheet. If False, 
    will use resdf as input. 


    Script will output images depending on output_tp argument. If 
    set to resamp, an image indicating the resampled pvalue at 
    each parcel will be generated. If output_tp set to cross, 
    images will be generated representing # of top hits, # of 
    appearances in thresholded results, and average statistic 
    (i.e. r or p value) at each parcel, respectively. Passing both 
    will generate all images. 

    Cross most appropriate if cross-validation performed with 
    threshold. CI most appropriate if jackknife or cross-validation 
    performed w/o threshold. resamp most appropriate if permute or 
    bootstrap performed w/o threshold.

    '''
    

    if par:
        resdf = parallel_in('coi','',par,'')
    
    # determine analysis resolutions
    sclz = []
    for ind in resdf.index.tolist():
        scl,jnk = ind.split(': ')
        resdf.ix[ind,'scale'] = scl
        resdf.ix[ind,'parcel'] = int(jnk)
        if scl not in sclz:
            sclz.append(scl)

    for scl in sclz:
        # prep inputs
        indf = resdf.loc[resdf['scale'] == scl]
        indf.index = indf[:]['parcel'].tolist()
        nind = indf.index.tolist()
        for i in range(int(scl)):
            if i not in indf.index.tolist():
                nind.append(i)
        indf = indf.reindex(nind,fill_value=0)
        for ind in indf.index.tolist():
            if not pandas.notnull(indf.ix[ind,'top_hits']):
                indf.ix[ind,'top_hits'] = 0
        indf = indf.sort()
        #### The line below should be improved.... ######
        scale_templ = glob(os.path.join(scldir,'%s%s.*'%(sclstr,scl)))[0]
        if output_tp == 'cross' or output_tp == 'all':
            print 'making top_hits image for scale %s'%(scl)
            wr.make_parcelwise_map(outdir,indf,scale_templ,
                    outfl=os.path.join(outdir,'%s_tophits_%s%s'%(outstr,input_tp,scl)),
                    add=True,col='top_hits')
            print 'making appearances image for scale %s'%(scl)
            wr.make_parcelwise_map(outdir,indf,scale_templ,
                    outfl=os.path.join(outdir,'%s_appearances_%s%s'%(outstr,input_tp,scl)),
                    add=True,col='appearances')
        if output_tp == 'ci' or output_tp == 'all':
                print 'making average coefficient image for scale %s'%(scl)
                wr.make_parcelwise_map(outdir,indf,scale_templ,
                        outfl=os.path.join(outdir,'%s_average_coef_%s%s'%(outstr,input_tp,scl)),
                        add=True,col='mean_%s'%(input_tp))
        if output_tp == 'resamp' or output_tp == 'all':
            if 'resample_pvalue' in indf.columns.tolist():
                print 'making resample pvalue image for scale %s'%(scl)
                wr.make_parcelwise_map(outdir,indf,scale_templ,
                        outfl=os.path.join(outdir,'%s_resample_p_%s%s'%(outstr,input_tp,scl)),
                        add=True,col='resample_pvalue')
            else:
                raise ValueError('resamp passed for argument output_tp, but no resampling statistics available')



