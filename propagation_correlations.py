import os
import nibabel as ni
import numpy as np
import pandas
import subprocess
import random,string
import scipy.stats.stats as st
import matplotlib.pyplot as plt
import jake_niak_utils as jni
from scipy.io import loadmat
from glob import glob

#image1 = '/Users/jakevogel/bellec_lab/aaic/density/masked_HCP_mean.nii.gz'
#image2 = '/Users/jakevogel/bellec_lab/aaic/density/spmT_0001.nii'
#msk = '/Users/jakevogel/bellec_lab/aaic/density/mask_GM.nii'
#xlab =
#ylab = 
#title =


def load_data(image1,image2,msk):
    mdata = ni.load(msk).get_data().flatten()
    # sometimes ni reads binary 1 as almost 1 floats so....
    dmask = []
    for x in mdata:
        if x != 0:
            dmask.append(1)
        else:
            dmask.append(0)

    dmask = np.logical_not(dmask)
    data1 = np.ma.masked_array(ni.load(image1).get_data(),mask = dmask).flatten()
    data2 = np.ma.masked_array(ni.load(image2).get_data(),mask = dmask).flatten()

    return data1,data2,dmask

def corr_n_plot(data1,data2):
    r,p = st.pearsonr(data1,data2)
    m,b = np.polyfit(data1,data2,1)
    plt.plot(data1,data2,'.')
    plt.plot(data1,m*data1+b,'-')
    print 'r = %s, p = %s'%(r,p)
    plt.show()

def codegen(N):
    cde = ''.join(random.choice(string.ascii_lowercase) for _ in range(N))

    return cde

def pthswp(wpth):
    oldpth = os.getcwd()
    os.chdir(wpth)

    return oldpth

def make_parcelwise_map(wpth,rdf,scale_templ, outfl, add = False, col = ''):
    oldpth = pthswp(wpth)
    cde = codegen(6)

    if not col:
        col = rdf.columns[0]
    for x in rdf.index.tolist():
        if add:
            i = x+1
        else:
            i = x
        val = rdf.ix[x,col]
        os.system('fslmaths %s -thr %s %s1'%(scale_templ,i,cde))
        os.system('fslmaths %s1.nii.gz -uthr %s %s2'%(cde,i,cde))
        os.system('fslmaths %s2.nii.gz -div %s2.nii.gz %s3'%(cde,cde,cde))
        if x == rdf.index.tolist()[0]:
            os.system('fslmaths %s3.nii.gz -mul %s %s'%(cde,val,outfl))
        else:
            os.system('fslmaths %s3.nii.gz -mul %s %s4'%(cde,val,cde))
            os.system('fslmaths %s4.nii.gz -add %s.nii.gz %s'%(cde,outfl,outfl))
        print 'seed %s built'%(i)
    os.system('rm %s*'%(cde))
    os.chdir(oldpth)

def img_2_maskd_array(img):
    msk = np.array(ni.load(img).get_data()).flatten()
    dmask = []
    for k in msk:
        if k != 0:
            dmask.append(1)
        else:
            dmask.append(0)

    msk = np.logical_not(dmask)

    return msk

def in_out_mask_ttest_from_conn_map(wpth,rseedz, imsk, omsk):

    oldpth = pthswp(wpth)

    incde = codegen(6)
    outcde = codegen(6)

    print 'generating data masks....'
    inmsk = img_2_maskd_array(imsk)
    outmsk = img_2_maskd_array(omsk)

    print 'beginning tests'
    ind = []
    for i in range(len(rseedz)):
        ind.append((i+1))
    df = pandas.DataFrame(np.zeros((len(rseedz),4)),index = ind, columns =['t','p','wt','wp'])
    for i,seed in enumerate(rseedz):
        os.system('fslmaths %s -mas %s %s'%(seed,imsk,incde))
        os.system('fslmaths %s -mas %s %s'%(seed,omsk,outcde))
        invals = np.ma.masked_array(ni.load('%s.nii.gz'%(incde)).get_data(),mask = inmsk).flatten()
        outvals = np.ma.masked_array(ni.load('%s.nii.gz'%(outcde)).get_data(),mask = outmsk).flatten()
        t,p = st.ttest_ind(invals,outvals)
        wt,wp = st.ttest_ind(invals,outvals,equal_var = False)
        df.ix[(i+1),'t'] = t 
        df.ix[(i+1),'p'] = p
        df.ix[(i+1),'wt'] = wt
        df.ix[(i+1),'wp'] = wp
        print 'finished with seed %s'%(seed)

    os.system('rm %s* %s*'%(incde,outcde))
    os.chdir(oldpth)

    return df

def in_out_mask_ttest_from_glm_file(wpth,scale_path,contrast,scale,imsk,omsk,parcel_img,membership=[],conndf=''):

    oldpth = pthswp(wpth)
    tdf = pandas.DataFrame(np.zeros((scale,4)), columns=['t','p','wt','wp'])

    if not membership:
        #determine which seeds are in which mask
        inseedz = []
        outseedz = []
        cde = codegen(6)
        print 'determining parcel membership...'

        for i in range(1,(scale+1)):
            os.system('fslmaths %s -thr %s -uthr %s %s'%(parcel_img,i,i,cde))
            os.system('fslmaths %s.nii.gz -mas %s %s1'%(cde,imsk,cde))
            os.system('fslmaths %s.nii.gz -mas %s %s2'%(cde,omsk,cde))
            ival = subprocess.check_output('fslstats %s1.nii.gz -V'%(cde),shell = True).rsplit()[0]
            oval = subprocess.check_output('fslstats %s2.nii.gz -V'%(cde),shell = True).rsplit()[0]
            # determine membership via winner-takes-all
            if int(ival) > int(oval):
                inseedz.append(i)
                print 'seed %s going inside mask'%(i)
            elif int(ival) < int(oval):
                outseedz.append(i)
                print 'seed %s going outside mask'%(i)
            else:
                print 'could not resolve seed %s. In vox = %s, out vox = %s. Excluding from analysis'%(i,ival,oval)
        os.system('rm %s*')
    else:
        inseedz = membership[0]
        outseedz = membership[1]


    print 'preparing connectivity map...'

    if type(conndf) == pandas.core.frame.DataFrame:
        df = conndf
    else:
        df = jni.create_df_from_glm(scale_path,contrast,pval=0.1)

    for i in range(scale):
        print 'calculating values for seed %s'%(i+1)
        ivalz = []
        ovalz = []
        indz = []
        for ind in df.index.tolist():
            for x in ind:
                if x == (i+1):
                    indz.append(ind)
        indz.remove(indz[i])
        for y in indz:
            if y[0] == (i+1):
                conn = y[1]
            else:
                conn = y[0]
            if conn in inseedz:
                ivalz.append(df.ix[y,'eff'])
            elif conn in outseedz:
                ovalz.append(df.ix[y,'eff'])
        invec = np.array(ivalz)
        outvec = np.array(ovalz)
        t,p = st.ttest_ind(invec,outvec)
        wt,wp = st.ttest_ind(invec,outvec,equal_var = False)
        tdf.ix[(i+1),'t'] = t 
        tdf.ix[(i+1),'p'] = p
        tdf.ix[(i+1),'wt'] = wt
        tdf.ix[(i+1),'wp'] = wp

    os.chdir(oldpth)

    return tdf,df,inseedz,outseedz

def get_parcelwise_correlation_map_from_connectivity_map(wpth,rseedz,parcel_img,nparcels):
    oldpth = pthswp(wpth)
    cde = codegen(6)
    df = pandas.DataFrame(np.zeros((len(nparcel),4)),columns = ['rval','rp','rho','rhop'])
    for num,seed in enumerate(rseedz):
        they = []
        for i in range(nparcels):
            os.system('fslmaths %s -thr %s -uthr %s %s'%(parcel_img,(i+1),(i+1),cde))
            os.system('fslmaths %s -mas %s.nii.gz %s1'%(cde,seed,cde))
            val1 = subprocess.check_output('fslstats %s1.nii.gz -M'%(cde), shell = True)
            they.append(val1)
        nthey = [float(x) for x in they]
        nthey = np.array(nthey)
        r,p = st.pearsonr(nthey,new_x)
        rho,rp = st.spearmanr(nthey,new_x)
        df.ix[num,'rval'] = r
        df.ix[num,'rp'] = p
        df.ix[num,'rho'] = rho
        df.ix[num,'rhop'] = rp
    os.system('rm %s*'%(cde))
    os.chdir(oldpth)

    return df

def get_parcelwise_correlation_map_from_glm_file(wpth,mat_path,scale_templ,scale,contrast,eff='r',cov_msk='',cov_img='',conndf= '',poly=1):
    if poly < 2 or type(poly) != int:
        print 'no polynomial contrast applied'
        poly = False

    oldpth = pthswp(wpth)
    rdf = pandas.DataFrame(np.zeros((scale,4)),columns =['rval','rp','rho','rhop'])

    if cov_msk != '':
        scalar = cov_msk
    else:
        if cov_img == '':
            raise IOError('you need to include either a covariate image or amaskd_np array')
        else:
            print 'resampling image covariate to parcels...'
            scalar = convert_voxels_to_parcels(wpth,cov_img,scale_templ,scale)
            scalar = [float(x) for x in scalar]
            scalar = np.array(scalar)

    if conndf == pandas.core.frame.DataFrame:
        df = conndf
    else:
        print 'converting matrix...'
        if eff == 'r':
            df = jni.create_df_from_mat(mat_path,scale,pval=0.1,mat_tp='glm')
            effcol = 'eff'
        elif eff == 't':
            df = jni.create_df_from_mat(mat_path,scale,pval=0.01,mat_tp='glm',eff_tp ='ttest')
            effcol = 'ttest'
        else:
            raise ValueError('a value of r (correlation map) or t (tmap) must be entered for eff')
    for i in range(scale):
        print 'calculating parcelwise correlation for seed %s'%(i)
        they = []
        for ind in df.index.tolist():
            for x in ind:
                if x == (i+1):
                    they.append(df.ix[ind,effcol])
        they.remove(they[i]) # removes a redundant connection
        nthey = np.array(they)
        r,p = st.pearsonr(nthey,scalar)
        rho,rp = st.spearmanr(nthey,scalar)
        if poly:
            poly_r1,pp1 = st.pearsonr(nthey,scalar**poly)
            poly_r2,pp2 = st.pearsonr(nthey**poly,scalar)
            poly_rho1,rpp1 = st.spearmanr(nthey,scalar**poly)
            poly_rho2,rpp2 = st.spearmanr(nthey**poly,scalar)
        rdf.ix[i,'rval'] = r
        rdf.ix[i,'rp'] = p
        rdf.ix[i,'rho'] = rho
        rdf.ix[i,'rhop'] = rp
        if poly:
            rdf.ix[i,'poly_rval1'] = poly_r1
            rdf.ix[i,'poly_rp1'] = pp1
            rdf.ix[i,'poly_rval2'] = poly_r2
            rdf.ix[i,'poly_rp2'] = pp2
            rdf.ix[i,'poly_rho1'] = poly_rho1
            rdf.ix[i,'poly_rhop1'] = rpp1
            rdf.ix[i,'poly_rho2'] = poly_rho2
            rdf.ix[i,'poly_rhop2'] = rpp2

    print 'finished all seeds'

    os.chdir(oldpth)

    return df, rdf, scalar

#def get_in_out_mask_vox_counts(percs, mskd_conmap,
#    for perc in percs:
#        os.system('fslmaths -thr %s clus'%(mskd_conmap,float(sorted(yconn)[perc])))
#        os.system('fslmaths clus.nii.gz -bin bin_clus')
#        tvox = subprocess.check_output('fslstats bin_clus.nii.gz -V',shell = True).rsplit()[0]
#        os.system('fslmaths sig0p05_spmt.nii.gz -mas bin_clus.nii.gz in_clus')
#        ivox = subprocess.check_output('fslstats in_clus.nii.gz -V',shell = True).rsplit()[0]
#        mvox = int(pvox) - int(ivox)
#        pdf.ix[perc,'tvox'] = int(tvox)
#        pdf.ix[perc,'ivox'] = int(ivox)
#        pdf.ix[perc,'mvox'] = int(mvox)

def convert_voxels_to_parcels(wpth,to_convert,scale_templ,scale):

    oldpth = pthswp(wpth)
    cde = codegen(6)
    scalar = []
    for i in range(scale):
        os.system('fslmaths %s -thr %s -uthr %s %s'%(scale_templ,(i+1),(i+1),cde))
        os.system('fslmaths %s -mas %s.nii.gz %s1'%(to_convert,cde,cde))
        val1 = subprocess.check_output('fslstats %s1.nii.gz -M'%(cde), shell = True)
        scalar.append(val1)

    os.system('rm %s*'%(cde))
    os.chdir(oldpth)

    return scalar

def generate_sliding_window_panel_from_matfiles(mat_pth,subdict,num_windows,scale,contrast,seed_connz,conconnz):
    """mat_pth = string path to where your connectivity matfiles live.
    subdict = a dictionary where subjects you want to include in the analysis are
    keys and the "window" they belong in is the value
    num_windows = int, the number of windows
    scale = int, what resolution the analysis was done in (i.e. how many parcels)
    contrast = string, what contrast to call from in the matfile (where the connectivity
    maps are)
    seed_conz = a list of labels for which you want to know the average
    connectivity
    conconnz = a list of labels representing for which parcels you want to find the
    average z(r) in relation to all parcells in seed_conz
    """

    items = []
    for i in range(num_windows):
        items.append('gp%s'%(i))

    pan = pandas.Panel(items=items,major_axis = subdict.keys(),minor_axis = conconnz)

    missing = []
    for sub in subdict.keys():
        conn = glob(os.path.join(mat_pth,'*%s*.mat'%(sub)))
        if len(conn) > 0:
            conn = conn[0]
            mat = loadmat(conn)[contrast][0][0][0][0]
            dicteff = jni.reshuffle_matrix(mat,scale)
            condf = pandas.DataFrame(np.zeros((len(mat),1)), index = sorted(dicteff.keys()), columns = ['z(r)'])
            for conx in condf.index.tolist():
                condf.ix[conx,'z(r)'] = dicteff[conx]
            for i in conconnz:
                valz = []
                for seed in seed_connz:
                    val = condf.ix[(i,seed),'z(r)']
                    if type(val) != np.float64:
                        val = condf.ix[(seed,i),'z(r)']
                    valz.append(val)
                avg = np.array(valz).mean()
                if type(avg) == np.float64:
                    pan.ix['gp%s'%(subdict[sub]),sub,i] = avg
                elif type(avg) == np.ndarray:
                    pan.ix['gp%s'%(subdict[sub]),sub,i] = avg[0]
            print 'just finished sub %s'%(sub)
        else:
            print 'could not find subject %s'%(sub)
            missing.append(sub)

    if missing:
        print 'the following subjects were not fount: %s'%(missing)

    return pan

def create_sliding_window_images(pan,scale_templ,outdir,outfl,add=False):

    for i in range(len(pan.items.tolist())):
        print 'working on window %s'%(i)
        ndf = pan.ix['gp%s'%(i)]
        connex = ndf.mean(axis=0).to_frame()
        make_parcelwise_map(outdir,connex,scale_templ,outfl=os.path.join(outdir,'%s%s'%(outfl,i)),add=add)

def identify_labels_within_mask(wpth,networks,mask,as_int=True):
    """networks is a list of paths to labeled atlases
    mask is a path to an image you wish to mask atlas with
    leave as_int true if you want label values to be returned as integers"""

    oldpth = pthswp(wpth)
    cde = codegen(6)

    labels_out = {}

    for net in networks:
        print 'working on network %s'%(net)
        os.system('fslmaths %s -mas %s %s'%(net,mask,cde))
        nimg = '%s.nii.gz'%(cde)
        data = np.array(ni.load(nimg).get_data().flatten())
        msk = np.logical_not(img_2_maskd_array(nimg))
        unique = []
        for vox in data[msk]:
            if vox not in unique:
                if as_int:
                    unique.append(int(vox))
                else:
                    unique.append(vox)

        labels_out.update({net: unique})

    os.system('rm %s'%(nimg))
    os.chdir(oldpth)

    return labels_out

def extract_data_from_specific_connections(labels,glm_key,outfl):
    """labels is a dict produced from the above function, with key = path to
    atlas and value = labels to extract
    glm_key is a dict mapping path to atlas (key) to a list where the first
    item is NIAK glm_file and the second is the number of parcels (value)
    outfl is label for output spreadsheet"""

    for atl,glm in glm_key.iteritems():
        scale_num = glm[1]
        glm = glm[0]
        print 'working on %s'%(glm)
        mat = loadmat(glm)
        p = mat['pce'][0]
        eff = mat['eff'][0]
        dictp = jni.reshuffle_matrix(p,int(scale_num))
        df = pandas.DataFrame(np.zeros((len(p),1)), index=sorted(dictp.keys()), columns = ['p'])
        for conn in df.index.tolist():
            df.ix[conn,'p'] = dictp[conn]

        targets = labels[atl]
        connz = []
        for lab in targets:
            for i in df.index.tolist():
                if lab in i:
                    connz.append(i)

        for conn in df.index.tolist():
            if conn not in connz:
                df.drop(conn,axix=0,inplace=True)

        flnm = atl.rsplit('.')[0]
        df.to_csv('%s_%s.csv'%(outfl,flnm))

def generate_rvals_maps_from_glm(wpth,scl_num,scale_templ,contr,outfl):

    oldpth = pthswp(wpth)
    cde = codegen(6)
    cnx = sorted(glob(os.path.join(wpth,'connectome*.mat')))

    for x in range(1,scl_num+1):
        os.system('fslmaths %s -thr %s -uthr %s -bin %s'%(scale_templ,x,x,'%s_%s'%(cde,x)))

    for cnt in cnx:
        pth,matnm = os.path.split(cnt)
        sid = matnm.rsplit('_')[1]
        print 'working on subject %s'%(sid)
        ofl = '%s_%s'%(sid,outfl)
        df = jni.create_df_from_mat(cnt,scl_no=scl_num,pval=0.05,mat_tp='ind',eff_tp=contr)
        d_seedz = {}
        for i in range(1,scl_num+1):
            seedz = []
            for conn in df.index.tolist():
                for c in conn:
                    if c == i:
                        if conn not in seedz:
                            seedz.append(conn)
            d_seedz.update({i:seedz})

        for seed,edges in d_seedz.iteritems():
            for edge in edges:
                if edge[0] == seed:
                    img_no = edge[1]
                else:
                    img_no = edge[0]
                val = df.ix[edge,df.columns[0]]
                tofl = '%s_in_%s'%(cde,img_no)
                os.system('fslmaths %s_%s.nii.gz -mul %s %s'%(cde,img_no,val,tofl))
            tmp = sorted(glob('*_in_*'))
            stg = 'fslmaths %s -add '%(tmp[0])
            for z in range(1,len(tmp)):
                if z == (len(tmp)-1):
                    stg = stg+'%s %s'%(tmp[z],'%s_3d_%s'%(cde,seed))
                else:
                    stg = stg+'%s -add '%(tmp[z])
            print 'building seedmap %s'%(seed)
            os.system(stg)

        volz = sorted(glob('*_3d_*'))
        vstg = 'fslmerge -t %s_scl%s_rmap'%(sid,scl_num)
        for y in range(len(volz)):
            vstg = vstg+' %s'%(volz[y])
        print 'building 4D volume...'
        os.system(vstg)

    os.system('rm %s*'%(cde))
    os.chdir(oldpth)
