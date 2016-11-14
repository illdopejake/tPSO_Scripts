import os,sys
from glob import glob
import pandas
import numpy as np
import statsmodels.formula.api as smf
from statsmodels.sandbox.stats.multicomp import multipletests as fwe
from statsmodels.sandbox.stats.multicomp import fdrcorrection0 as fdr
sys.path.insert(0,'/Users/jakevogel/git/tPSO_scripts/')
import propagation_correlations as wr

srch_str = 's*.txt'
varsheet = '/Users/jakevogel/Desktop/tPSO.csv'

# Make sure the first variable is your independent variable of interest.
# Also make sure each string corresponds directly to a column in varsheet
varlist = ['tPSO','Age','Gender','FD','ApoE4_status']

def parcelwise_analysis(dirz,srch_str,varsheet,varlist):
    res = {}
    sigz = {}
    for direc in dirz:
        jnk,ref = os.path.split(direc)
        print 'working on directory %s'%(ref)
        print 'forming matrix'
        mapz = sorted(glob(os.path.join(direc,srch_str)))
        jnk = pandas.read_table(mapz[0],header=None)
        jnk.drop(jnk.columns[-1],axis=1,inplace=True)
        sz = len(jnk)
        mtx = np.full((sz,len(mapz)),np.nan)
        for i,map in enumerate(mapz):
            mdf = pandas.read_table(map,header=None)
            mdf.drop(mdf.columns[-1],axis=1,inplace=True)
            denz = mdf.sum(axis=1)
            for x,d in enumerate(denz.tolist()):
                mtx[x,i] = d

        print 'creating spreadsheet'
        if varsheet.split('.')[1][-3:] == 'csv':
            cdf = pandas.read_csv(varsheet)
        else:
            cdf = pandas.ExcelFile(varsheet).parse('Sheet1')

        if len(cdf) != mtx.shape[1]:
            raise IOError('input varsheet must have the same number of rows as there are text files in your directories')
        for v in varlist:
            if v not in cdf.columns.tolist():
                raise IOError('all items in varlist must correspond to columns in varsheet')

        for j,sub in enumerate(cdf.index.tolist()):
            for i in range(len(mtx)):
                cdf.ix[sub,'p%s_dens'%(i)] = mtx[i,j]

        print 'running models'
        rdf = pandas.DataFrame(np.full((len(mtx),2),np.nan),
            columns =['t','p'])

        for i in range(len(mtx)):
            stmnt = build_statement('p%s_dens'%(i),varlist)
            lm = smf.ols(stmnt,data=cdf).fit()
            rdf.ix[i+1,'t'] = lm.tvalues[1]
            rdf.ix[i+1,'p'] = lm.pvalues[1]

        rdf.drop(rdf.index[0],axis=0,inplace=True)

        print 'correcting models'
        fdrtst = fdr(np.array(rdf[:]['p'].tolist()))
        fwetst = fwe(np.array(rdf[:]['p'].tolist()))
        for i in range(len(fdrtst[1])):
            rdf.ix[i+1,'fdr'] = fdrtst[1][i]
        for i in range(len(fwetst[1])):
            rdf.ix[i+1,'fwe'] = fwetst[1][i]

        res.update({ref: rdf})

        sig = []
        for parc in rdf.index.tolist():
            if rdf.ix[parc,'fdr'] < 0.1 or rdf.ix[parc,'fwe'] < 0.1:
                sig.append((parc,rdf.ix[parc,'fdr'],rdf.ix[parc,'fwe']))

        sigz.update({ref: sig})

    return res,sigz

def build_statement(dvar,varlist):
    stmnt = '%s ~ %s +'%(dvar,varlist[0])
    for i,v in enumerate(varlist):
        if i != 0:
            if i != (len(varlist) - 1):
                stmnt = '%s %s +'%(stmnt,v)
            else:
                stmnt = '%s %s'%(stmnt,v)
    return stmnt

def create_output_images(indict, scale_templ, cols, outdir):

    for ref,rdf in indict.iteritems():
        for col in cols:
            outfl = os.path.join(outdir,'parcelwise_%s_%smap'%(ref,col))
            wr.make_parcelwise_map(outdir,rdf,scale_templ,outfl,add=False,col=col)

def extract_network_densities(dirz,srch_str,subs,netsht):
    # assumes netdf is a spreadsheet with parcels as index and network
    # membership as the first columns

    print 'identifying networks'
    if netsht.split('.')[1][-3:] == 'csv':
        netdf = pandas.read_csv(netsht)
    else:
        netdf = pandas.ExcelFile(netsht).parse('Sheet1')
    netz = []
    for parc in netdf.index.tolist():
        if netdf.ix[parc,netdf.columns[0]] not in netz:
            if  pandas.notnull(netdf.ix[parc,netdf.columns[0]]):
                netz.append(netdf.ix[parc,netdf.columns[0]])
    print netz
    cdict = {}
    for net in netz:
        netno = []
        for parc in netdf.index.tolist():
            if netdf.ix[parc,netdf.columns[0]] == net:
                netno.append('1')
        cdict.update({net: len(netno)})

    if netdf.index.tolist()[0] == 1:
        nind = [x-1 for x in netdf.index.tolist()]
        netdf.index = nind

    den_dfz = {}
    for direc in dirz:
        jnk,ref = os.path.split(direc)
        print 'working on directory %s'%(ref)
        mapz = sorted(glob(os.path.join(direc,srch_str)))
        cdf = pandas.DataFrame(index=subs)
        for i,sub in enumerate(subs):
            print 'working on subject %s'%(sub)
            mdf = pandas.read_table(mapz[i],header=None)
            mdf.drop(mdf.columns[-1],axis=1,inplace=True)
            for net in netz:
                catch = []
                for x in mdf.index.tolist():
                    if netdf.ix[x,'membership'] == net:
                        for y in mdf.index.tolist():
                            if netdf.ix[y,'membership'] == net:
                                if mdf.ix[x,y] == 1:
                                    catch.append((x,y))
                netdens = float(len(catch)) / (cdict[net] * (cdict[net] - 1))
                cdf.ix[sub,'%s_density'%(net)] = netdens
        den_dfz.update({ref: cdf})

    return den_dfz

def directly_compare_coordinates(dirz,srch_str,subs,xlist,ylist,yname):
    den_dfz = {}
    for direc in dirz:
        jnk,ref = os.path.split(direc)
        print 'working on directory %s'%(ref)
        mapz = sorted(glob(os.path.join(direc,srch_str)))
        cdf = pandas.DataFrame(index=subs)
        for i,sub in enumerate(subs):
            print 'working on subject %s'%(sub)
            mdf = pandas.read_table(mapz[i],header=None)
            mdf.drop(mdf.columns[-1],axis=1,inplace=True)
            for x in mdf.index.tolist():
                catch = 0
                if x in xlist:
                    for y in mdf.index.tolist():
                        if y in ylist:
                            if mdf.ix[x,y] == 1:
                                catch = catch + 1
                    cdf.ix[sub,'p%s_to_%s'%(x,yname)] = catch
        den_dfz.update({ref: cdf})

    return den_dfz
