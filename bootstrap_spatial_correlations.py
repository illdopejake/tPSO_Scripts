import os
from glob import glob
import pandas
import random
import numpy as np
import nibabel as ni
import scipy.stats.stats as st
import propagation_correlations as wr

def define_bootstrap_sample(subdf,subcol,subpth,pv,sample_perc=0.5,num_gps=3,taskid = ''):
"""takes an existing sample of subjects from a spreadsheet and creates a
subsample, balanced by a variable. Outputs a subsamble membership txt file and
a list of paths corresponding to subsample member scans)

subdf = a spreadsheet containing at least a list of subjects and values for a
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

    if sample_perc > 1 or if sample_perc <=0:
        raise IOError('sample_perc must be a float between 0 and 1')

    if subdf[-3:] == 'xls' or subdf[-4:] == 'xlsx':
        subdf = pandas.ExcelFile(subdf).parse('Sheet1')
    elif subdf[-3:] == 'csv':
        subdf = pandas.read_csv(subdf)
    else:
        raise IOError('input spreadsheet filetype not recognized. Please use
        .xls, .xlsx or .csv')

    if subcol != 'index':
        if subcol in subdf.columns.tolist():
            df.index = df[:][subcol].tolist()
        else:
            raise IOError('there is no column in the spreadsheet called
            %s'%(subcol))

    for sub in subdf.index.tolist():
        if len(glob(os.path.join(subpth,'%s*'%(sub)))[0]) < 1:
            print 'no scan was found for subject %s. Removing from
            analysis.'%(sub)
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

   # collect scans 
    scans = []
    for sub in subsamp:
        scn = glob(os.path.join(subpth,'*%s*'%(sub)))[0]
        scans.append(scn)

    # make membership record
    ndf = pandas.DataFrame(index=subsamp)
    if len(taskid) > 0:
        ndf.to_csv('%s_subsample_membership'%(taskid))
    else:
        cde = wr.codegen(6
        ndf.to_csv('%s_subsample_membership'%(cde))

    return scans

def voxelwise_analysis

def spatial_correlation_searchlight

def id_sig_results

def replicate_previous_findings
