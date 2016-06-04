import os
from glob import glob
from scipy.io import loadmat
import pandas
import numpy as np

def reshuffle_matrix(array,scale_n,var=False):
    pdict = {}
    const = scale_n
    cntr = 1
    bcounter = 0

    for x in range(scale_n):
        while cntr < const+1:
            if var:
                val = array[bcounter]
                pdict.update({(x+1,cntr+x): val[0]})
            else:
                pdict.update({(x+1,cntr+x): array[bcounter]})
            cntr = cntr+1
            if bcounter != len(array):
                bcounter = bcounter + 1
            else:
                cntr = const+1
        const = const - 1
        cntr = 1

    return pdict

def create_df_from_glm(scale,contrast,pval,eff_tp='eff',addition=False,addno=0):

    pth,nme = os.path.split(scale)
    if nme == '':
        pth,nme = os.path.split(pth)
    jnk,scale_num = nme.split('ale')


    if addition:
        if int(scale_num) > 6:
            scale_num = int(scale_num)+int(addno)

    cont = os.path.join(scale,contrast)
    mat = loadmat(os.path.join(cont,'glm_%s_%s.mat'%(contrast,nme)))

    if eff_tp not in mat.keys():
        raise ValueError('please set eff_tp to a valid key of the glm.mat file')

    p = mat['pce'][0]
    eff = mat[eff_tp][0]

    dictp = reshuffle_matrix(p,int(scale_num))
    dicteff = reshuffle_matrix(eff,int(scale_num))

    df = pandas.DataFrame(np.zeros((len(p),5)), index = sorted(dictp.keys()), columns = [eff_tp,'p','fdr','sig','fdr_sig'])

    ddict = {eff_tp: dicteff, 'p': dictp}

    for lab,dick in ddict.iteritems():
        for conn in df.index.tolist():
            df.ix[conn,lab] = dick[conn]

    fdr = mat['fdr']
    for ind in df.index.tolist():
        df.ix[ind,'fdr'] = fdr[(ind[0]-1)][(ind[1]-1)]

    fdrs = mat['test_q']
    for ind in df.index.tolist():
        df.ix[ind,'fdr_sig'] = fdrs[(ind[0]-1)][(ind[1]-1)]

    for k,v in dictp.iteritems():
        if v < pval:
            df.ix[k, 'sig'] = 1

    return df

def save_sig_results(df,corr='fdr',eff_tp = 'eff'):
    
    if corr == 'fdr':
        ccol = 'fdr_sig'
    else:
        ccol = 'sig'

    inds = df.index.tolist()
    ninds = [str(x) for x in inds]
    df.index = ninds

    sig_cols = []

    for sub in df.index.tolist():
        if df.ix[sub,ccol] == 1:
            sig_cols.append(sub)

    sigdf = pandas.DataFrame(np.zeros((len(sig_cols),2)),index = sig_cols, columns = [eff_tp,'p'])
    
    for sub in df.index.tolist():
        if df.ix[sub,ccol] == 1:
            sigdf.ix[sub,eff_tp] = df.ix[sub][1]
            if corr == 'fdr':
                sigdf.ix[sub,'p'] = df.ix[sub,'fdr']
            else:
                sigdf.ix[sub,'p'] = df.ix[sub,'p']

    return sigdf

def determine_top_connections(df,seed,typ = 'fdr',perc=0.01):
    
    if typ not in df.columns.tolist():
        raise ValueError('Please set typ to a valid column name, e.g. fdr or eff') 

    allconz = []

    for conn in df.index.tolist():
        if seed in conn:
            allconz.append(df.ix[conn,typ])

    scale = len(allconz)
    est_seed = float(scale) * perc
    if (int(est_seed) + .5) > est_seed:
        no_seeds = int(est_seed) 
    else:
        no_seeds = int(est_seed) + 1

    if typ == 'eff' or typ == 'std_eff' or typ == 'ttest':
        top_ps = sorted(allconz)[-(no_seeds):]
    else:
        top_ps = sorted(allconz)[:no_seeds]
    
    top_dict = {}

    for conn in df.index.tolist():
        if seed in conn:
            if df.ix[conn,typ] in top_ps:
                top_dict.update({conn: df.ix[conn,typ]})
    
    top_cs = pandas.DataFrame(top_dict,index = ['p'])

    return top_cs   
